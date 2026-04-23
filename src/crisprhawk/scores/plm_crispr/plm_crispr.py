"""PLM-CRISPR on-target efficiency scoring for CRISPR-HAWK.

This module provides a self-contained inference interface for the PLM-CRISPR
model, which predicts on-target CRISPR editing efficiency by jointly encoding
the sgRNA sequence and the Cas9 protein structure via ESM per-residue
embeddings.

The module is intentionally stateless: the model and protein embeddings are
loaded fresh on each call to ``compute_plm_crispr_score``, keeping the
function safe for use in multiprocessing contexts
"""

from torch import device, Tensor
from typing import List

import torch.nn as nn

import importlib
import torch
import warnings
import sys
import os

# pass BATCH_SIZE sequences at a time to avoid oom errors
BATCH_SIZE = 64

# define fixed spcas9 scaffold sequence (first 36 nt), used to extend the 23-nt
# spcer+PAM to the 59-nt input expected by the PLM-CRISPR model
_SCAFFOLD = (
    "gtttCagagctaTGCTGgaaaCAGCAtagcaagttGaaataaggctagtcc"
    "gttatcaacttgaaaaagtggcaccgagtcggtgcttt"
)[:36].upper()

# define expected total input length for PLM-CRISPR (spacer+PAM + scaffold)
_SEQLEN = 59

# nucleotide - integer encoding used by PLM-CRISPR's data loader
_TRANSDIR = {nt: i for i, nt in enumerate(["A", "T", "C", "G", "N"])}

# define mapping from CRISPR-HAWK guessed cas system to the corresponding PLM-CRISPR
# protein embedding filename and dataset-type index.
# SPCAS9 = 3, XCAS9 = 4 (only systems scored by PLM-CRISPR in CRISPR-HAWK).
# typ index follows the dataset ordering used during PLM-CRISPR training:
#   6 -> WT_kim / WT_Wang (SpCas9 WT), 10 -> XCas9
_CAS_CONFIG = {
    3: {
        "protein_file": "SpCas9.esm.pt",
    },  # SPCAS9
    4: {
        "protein_file": "xCas9.esm.pt",
    },  # XCAS9
}

# define PLM-CRISPR models folder
_MODELSDIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "models")


def _load_model(device: device) -> nn.Module:
    """Load the PLM-CRISPR model from the original full-object checkpoint.

    The checkpoint ``100.pt`` was serialized with ``torch.save(net, ...)``
    in the original PLM-CRISPR training environment, where ``train_general``
    and ``Protein`` were importable as top-level modules. To deserialize it
    correctly in CRISPR-HAWK — where these modules live inside the
    ``crisprhawk.scores.plm_crispr`` package — this function temporarily
    registers them under their bare names in ``sys.modules`` so that pickle's
    ``find_class`` can resolve the stored class references.

    The ``sys.modules`` dictionary is always restored to its original state
    after loading, whether or not an exception occurs, to avoid polluting the
    global import namespace.

    Args:
        device (torch.device): Device to map the loaded model onto (CPU or
            GPU).

    Returns:
        nn.Module: The PLM-CRISPR model in eval mode, ready for inference.

    Raises:
        ImportError: If ``train_general`` or ``Protein`` cannot be imported
            as part of the ``crisprhawk.scores.plm_crispr`` package.
    """
    # PLM-CRISPR models are stored as full serialized objects (not state_dicts)
    model_path = os.path.join(_MODELSDIR, "100.pt")
    # Import train_general and Protein as proper package modules so that their
    # internal relative imports (e.g. from .Protein import TextCNN) are
    # resolved correctly, then register them under the bare names that pickle
    # stored in the checkpoint
    _pkg = "crisprhawk.scores.plm_crispr"
    _shims = {
        "train_general": importlib.import_module(f"{_pkg}.train_general"),
        "Protein": importlib.import_module(f"{_pkg}.Protein"),
    }
    # Snapshot the current state of sys.modules for the two keys so we can
    # restore it exactly after loading (None means the key was absent)
    _prev = {k: sys.modules.get(k) for k in _shims}
    sys.modules.update(_shims)
    try:
        net = torch.load(model_path, map_location=device, weights_only=False)
    finally:
        # restore sys.modules uncoditionally to avoid import namespace leaks
        for k, v in _prev.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    net = net.to(device)
    net.eval()
    return net


def _load_esm_embedding(protein_file: str, device: device) -> Tensor:
    """Load the ESM per-residue protein embedding for a given Cas9 variant.

    Args:
        protein_file (str): Filename of the protein embedding.
        device (torch.device): Device to load the tensor onto.

    Returns:
        torch.Tensor: Per-residue ESM embedding of shape ``(L, 1280)``,
            cast to float16 to match the dtype expected by the PLM-CRISPR
            protein branch.

    Raises:
        KeyError: If the loaded checkpoint does not contain a
            ``representations`` key.
    """
    protein_file_path = os.path.join(_MODELSDIR, protein_file)
    protein_data = torch.load(
        protein_file_path, map_location=device, weights_only=False
    )
    return protein_data["representations"].to(torch.float16)


def _concat_scaffold(grna: str) -> str:
    """Append the fixed SpCas9 scaffold to a spacer+PAM string.

    Concatenates the guide's spacer+PAM with the 36-nt scaffold and
    uppercases the result. Raises if the resulting length is not exactly
    59 nt, which would indicate an unexpected guide length.

    Args:
        grna (str): Spacer+PAM sequence (expected 23 nt for SpCas9).

    Returns:
        str: Uppercase 59-nt sequence ready for encoding.

    Raises:
        ValueError: If the concatenated sequence length differs from 59.
    """
    s = (grna + _SCAFFOLD).upper()
    if len(s) != _SEQLEN:
        raise ValueError(
            f"PLM-CRISPR input sequence has length {len(s)}, expected "
            f"{_SEQLEN}. Guide='{grna}' (len={len(grna)})"
        )
    return s


def _build_input_sequences(guides: List[str]) -> List[str]:
    """Builds 59-nt PLM-CRISPR input sequences from Guide objects.

    Concatenates each guide's spacer+PAM sequence (guidepam, 23 nt for SpCas9)
    with the fixed 36-nt scaffold to reach the required 59-nt input length.
    Sequences are uppercased before encoding.

    Args:
        guides (List[str]): List of guides. Each guide contains (PAM+spacer or
            spacer+PAM depending on orientation).

    Returns:
        List[str]: List of 59-nt uppercase strings ready for encoding.

    Raises:
        ValueError: If any resulting sequence is not exactly 59 nt.
    """
    return [_concat_scaffold(g) for g in guides]


def _encode_sequences(sequences: List[str]) -> Tensor:
    """Encode a list of 59-nt nucleotide strings into a LongTensor.

    Maps each nucleotide to its integer index using the PLM-CRISPR encoding:
    A=0, T=1, C=2, G=3, N=4. The resulting tensor has shape ``(N, 59)``
    where N is the number of sequences.

    Args:
        sequences (List[str]): Uppercase 59-nt strings to encode.

    Returns:
        torch.Tensor: LongTensor of shape ``(N, 59)`` containing integer-
            encoded nucleotide sequences.

    Raises:
        KeyError: If any character in a sequence is not present in
            ``_TRANSDIR`` (i.e. not one of A, T, C, G, N).
    """
    encoded = [[_TRANSDIR[nt] for nt in seq] for seq in sequences]
    return torch.LongTensor(encoded)


def compute_plm_crispr_score(guides: List[str], cas_system: int) -> List[float]:
    # retrieve cas-system-specific config
    config = _CAS_CONFIG[cas_system]
    # resolve device (GPU or CPU)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    net = _load_model(device)  # load plm model
    # load esm protein embedding for the input cas variant
    protein_emb = _load_esm_embedding(config["protein_file"], device)
    # build and encode input guides
    sequences = _build_input_sequences(guides)
    sgrna_tensor = _encode_sequences(sequences).to(device)
    B = sgrna_tensor.shape[0]  # retrieve batch size
    # expand protein embedding to match batch size
    # protein_batch = protein_emb.unsqueeze(0).repeat(B, 1, 1).to(device)
    # perform inference
    sgrna_tensor = sgrna_tensor.long().to(device)
    # use dynamic batch inference
    batch_size = BATCH_SIZE
    scores: List[float] = []  # plm-crispr scores
    while batch_size >= 1:  # initialize batch size to constant
        scores.clear()  # clear scores on retry to prevent appending duplicates
        try:
            for start_idx in range(0, B, batch_size):
                end_idx = min(B, start_idx + batch_size)
                current_batch_size = end_idx - start_idx
                # slice the current batch of guides
                sgrna_batch = sgrna_tensor[start_idx:end_idx, :]
                # expand protein embeddings to match the current batch size
                protein_batch = protein_emb.unsqueeze(0).repeat(current_batch_size, 1, 1).to(device)
                with torch.no_grad():
                    predictions = net(sgrna_batch, 0, protein_batch, train=False)
                # flatten predictions and append to scores
                scores_batch = predictions.view(-1).cpu().tolist()
                scores.extend(scores_batch)
            # if the loop completes without throwing an error, exit
            break
        except RuntimeError as e:
            if "out of memory" not in str(e).lower():
                # raise if the RuntimeError is not related to memory
                raise
            # free up VRAM before attempting the smaller batch size
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            new_batch_size = batch_size // 2
            # if we OOM on batch size 1, we can't go any lower
            if new_batch_size == 0:
                raise RuntimeError(
                    "CUDA Out of Memory even with batch size 1."
                    "Retry by setting CUDA_VISIBLE_DEVICES=-1 to run on CPU."
                ) from e
            warnings.warn(f"OOM Detected: Retrying with batch size = {new_batch_size}")
            batch_size = new_batch_size
    # ensure output is always a list even for a single guide
    return [scores] if isinstance(scores, float) else scores
