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
import sys
import os


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
    3: {"protein_file": "SpCas9.esm.pt", },   # SPCAS9
    4: {"protein_file": "xCas9.esm.pt",  },  # XCAS9
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
        "Protein":       importlib.import_module(f"{_pkg}.Protein"),
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

    ESM embeddings are stored as PyTorch checkpoint files containing a
    ``representations`` key with the per-residue embedding tensor of shape
    ``(L, 1280)``, where L is the protein sequence length and 1280 is the
    ESM-1b embedding dimension.

    Args:
        protein_file (str): Filename of the protein embedding (e.g.
            ``SpCas9.esm.pt``), looked up inside ``_MODELSDIR``.
        device (torch.device): Device to load the tensor onto.

    Returns:
        torch.Tensor: Per-residue ESM embedding of shape ``(L, 1280)``,
            cast to float16 to match the dtype expected by the PLM-CRISPR
            protein branch.

    Raises:
        FileNotFoundError: If the protein embedding file does not exist.
        KeyError: If the loaded checkpoint does not contain a
            ``representations`` key.
    """
    protein_file_path = os.path.join(_MODELSDIR, protein_file)
    protein_data = torch.load(protein_file_path, map_location=device, weights_only=False)
    return protein_data["representations"].to(torch.float16)

def _concat_scaffold(grna: str) -> str:
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
    """
    return [_concat_scaffold(g) for g in guides]

def _encode_sequences(sequences: List[str]) -> Tensor:
    """Encodes a list of 59-nt sequences into a LongTensor for PLM-CRISPR.

    Each nucleotide is mapped to an integer index using the PLM-CRISPR base
    encoding (A=0, T=1, C=2, G=3, N=4). The output tensor has shape
    (N, 59) where N is the number of sequences.

    Args:
        sequences (List[str]): List of 59-nt uppercase sequences.

    Returns:
        torch.Tensor: LongTensor of shape (N, 59) with encoded sequences.

    Raises:
        KeyError: If any nucleotide is not in the encoding dictionary.
    """
    
    for s in sequences:
        print(s)

    encoded = [[_TRANSDIR[nt] for nt in seq] for seq in sequences]
    return torch.LongTensor(encoded)


def compute_plm_crispr_score(guides: List[str], cas_system: int) -> List[float]:
    """Run PLM-CRISPR inference for a list of guides of the same Cas system.

    Loads the serialized PLM-CRISPR model and the corresponding ESM protein
    embedding for the given Cas system, builds the 59-nt input sequences,
    encodes them, and returns predicted on-target efficiency scores.

    The model and protein embedding are loaded fresh each call (no global
    state), which keeps the function stateless and safe for multiprocessing.
    GPU is used automatically if available, otherwise CPU.

    Args:
        guides (List[Guide]): List of Guide objects to score. All guides must
            belong to the same Cas system.
        cas_system (int): CRISPR-HAWK Cas system identifier (SPCAS9=3,
            XCAS9=4), used to look up the correct protein file and typ index.

    Returns:
        List[float]: Predicted PLM-CRISPR efficiency scores, one per guide,
            in the same order as the input list.

    Raises:
        KeyError: If ``cas_system`` is not in ``_PLMCRISPR_CAS_CONFIG``.
        FileNotFoundError: If the model or protein embedding file is missing.
        ValueError: If any input sequence cannot be constructed at 59 nt.
    """
    # retrieve cas-system-specific config
    config = _CAS_CONFIG[cas_system]
    
    typ_idx = 0
    # resolve device (GPU or CPU)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    net = _load_model(device)

    # load esm protein embedding for the input cas variant
    protein_emb = _load_esm_embedding(config["protein_file"], device)

   

    # # build and encode input guides
    sequences = _build_input_sequences(guides)
    sgrna_tensor = _encode_sequences(sequences).to(device)

    print(f"shape: {sgrna_tensor.shape}")
    print(f"{sgrna_tensor[0]}")

    batch_size = sgrna_tensor.shape[0]
    # batch_size = sgrna_tensor.shape[0]  # batch size

    # expand protein embedding to match batch size
    protein_batch = protein_emb.unsqueeze(0).repeat(batch_size, 1, 1).to(device)

    typ_tensor = typ_idx  # the model accepts a plain int for typ
    
    # inference
    sgrna_tensor = sgrna_tensor.long().to(device)
    with torch.no_grad():
        predictions = net(sgrna_tensor, typ_idx, protein_batch, train=False)
    
    # retrieve scores
    scores = predictions.squeeze().cpu().tolist()
    # ensure output is always a list even for a single guide
    if isinstance(scores, float):
        scores = [scores]
    return scores


