"""Compute sgDesigner on-target activity scores for CRISPR guides.

This module prepares guide sequences, invokes the external sgDesigner tool,
and parses its output to return per-guide scoring results.
"""

from ...utils import create_folder

from typing import List, Tuple, Dict

import os
import subprocess
import tempfile


def _find_output_txt(sgdesigner_outdir: str) -> str:
    """Locate the sgDesigner output text file within the result directory.

    This function searches recursively for the expected sgDesigner result file
    and returns the first matching path it finds.

    Args:
        sgdesigner_outdir: Path to the directory where sgDesigner wrote its output.

    Returns:
        The full path to the sgDesigner prediction result text file.
    """
    candidates = []
    for root, _, files in os.walk(sgdesigner_outdir):
        candidates.extend(
            os.path.join(root, fname)
            for fname in files
            if fname == "sgDesigner_V2.0_prediction_result.txt"
        )
    assert bool(candidates)  # sgdesigner failed?
    return candidates[0]


def _load_sgdesigner_scores(txt_path: str, expected_26mers: List[str]) -> List[float]:
    """Load sgDesigner scores from a text result file and align them to expected
    guides.

    This function parses sgDesigner output lines, matches each score to the
    corresponding submitted spacer sequence, and returns scores ordered by guide
    index.

    Args:
        txt_path: Path to the sgDesigner prediction result text file.
        expected_26mers: List of expected 26-nt guide sequences in submission order.

    Returns:
        A list of sgDesigner scores corresponding to each expected guide in
            input order.

    Raises:
        RuntimeError: If any expected guide is missing an sgDesigner score in the
            result file.
    """
    scores: List[float] = [float("nan")] * len(expected_26mers)
    expected_map = {
        f"guide_{i}": seq[:20].upper() for i, seq in enumerate(expected_26mers)
    }
    with open(txt_path, mode="r") as fin:
        fin.readline()  # skip header
        for line in fin:
            gid, score_raw, spacer, _, _ = line.split("\t")
            guide_id = gid.split("_")[1]
            if f"guide_{guide_id}" not in expected_map:
                continue
            # keep only the row matching the submitted 20-nt spacer
            if spacer.upper() != expected_map[f"guide_{guide_id}"]:
                continue
            scores[int(guide_id)] = float(score_raw)
    missing = [i for i, s in enumerate(scores) if s != s]  # NaN check
    if missing:
        raise RuntimeError(
            f"Missing sgDesigner scores for {len(missing)} guides. "
            f"First missing indices: {missing[:10]}"
        )
    return scores


def _generate_sgdesigner_tmp_data(tmpdir: str) -> Tuple[str, str, str]:
    """Create sgDesigner input and output paths inside a temporary directory.

    This function defines where the guide FASTA file, result files, and
    temporary working files will be stored for a single sgDesigner run.

    Args:
        tmpdir: Path to the temporary base directory where sgDesigner data
            will be created.

    Returns:
        A tuple containing the paths to the guides FASTA file, the sgDesigner
        results directory, and the temporary working directory.
    """
    # generate guides fasta required and out folder by crispron script
    return (
        os.path.join(tmpdir, "guides.fa"),
        os.path.join(tmpdir, "sgdesigner_results"),
        os.path.join(tmpdir, "temp"),
    )


def _write_guides_fasta(guides: List[str], sgdesigner_fasta: str) -> None:
    """Write guide sequences to a FASTA file in the format expected by sgDesigner.

    This function assigns a unique identifier to each guide and saves them
    in uppercase FASTA format for downstream sgDesigner processing.

    Args:
        guides: List of guide sequences to write to the FASTA file.
        sgdesigner_fasta: Path to the FASTA file where guides will be written.
    """
    with open(sgdesigner_fasta, mode="w") as fout:
        for i, seq in enumerate(guides):
            fout.write(f">guide_{i}\n{seq.upper()}\n")
    assert os.stat(sgdesigner_fasta).st_size > 0


def compute_sgdesigner_score(
    guides: List[str], conda: str, env_name: str
) -> List[float]:
    """Compute sgDesigner on-target activity scores for a list of guide sequences.

    This function runs the external sgDesigner tool in a specified conda
    environment and returns one score per input guide in the original order.

    Args:
        guides: List of guide sequences to be scored by sgDesigner.
        conda: Path to the conda executable used to invoke the sgDesigner environment.
        env_name: Name of the conda environment in which sgDesigner is installed.

    Returns:
        A list of sgDesigner scores corresponding to each guide in the input list.

    Raises:
        RuntimeError: If the sgDesigner script fails or does not produce scores
            for all provided guides.
    """
    assert bool(guides)  # otherwise we shouldn't be here
    # get path to sgdesigner script
    sgdesigner_root = os.path.abspath(os.path.dirname(__file__))
    sgdesigner_pl = os.path.join(sgdesigner_root, "sgDesigner.pl")
    assert os.path.isfile(sgdesigner_pl)
    # score guides with sgdesigner
    # with tempfile.TemporaryDirectory(prefix="sgdesigner_", dir=sgdesigner_root) as tmpdir:

    tmpdir = tempfile.mkdtemp(prefix="sgdesigner_")
    # generate guides fasta, output folder, and tmp folder required
    # by sgdesigner's script
    sgdesigner_fasta, sgdesigner_results, sgdesigner_tmpdir = (
        _generate_sgdesigner_tmp_data(tmpdir)
    )
    for d in [sgdesigner_results, sgdesigner_tmpdir]:
        create_folder(d, exist_ok=True)  # create sgdesigner folders
    # write guides to sgdesigner fasta
    _write_guides_fasta(guides, sgdesigner_fasta)
    # command to call sgDesigner perl script
    cmd = (
        f"export SGDESIGNER_RESULT_DIR='{sgdesigner_results}'; "
        f"export SGDESIGNER_TEMP_DIR='{sgdesigner_tmpdir}'; "
        f"perl '{sgdesigner_pl}' -f '{sgdesigner_fasta}'"
    )
    try:
        subprocess.run(
            [conda, "run", "-n", env_name, "bash", "-c", cmd],
            check=True,  # Raise exception on non-zero exit code
            capture_output=True,  # Capture stderr for better debugging
            text=True,
            cwd=sgdesigner_root,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"sgDesigner script failed with exit code {e.returncode}: {e.stderr}"
        ) from e
    txt_path = _find_output_txt(sgdesigner_results)
    return _load_sgdesigner_scores(txt_path, guides)
