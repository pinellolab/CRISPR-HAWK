"""Compute CRISPRon on-target activity scores for CRISPR guides.

This module prepares guide sequences, invokes the external CRISPRon tool,
and parses its output to return per-guide scoring results.
"""

from typing import List, Tuple

import csv
import os
import subprocess
import tempfile


def _find_output_csv(crispron_tmpdir: str) -> str:
    """Locate the CRISPRon output CSV file within a temporary directory.

    This function searches recursively for CSV files and prefers common
    CRISPRon result filenames when multiple candidates are found.

    Args:
        crispron_tmpdir: Path to the temporary directory produced by CRISPRon.

    Returns:
        The path to the selected CSV file containing CRISPRon results.
    """
    candidates = []
    for root, _, files in os.walk(crispron_tmpdir):
        candidates.extend(
            os.path.join(root, fname)
            for fname in files
            if fname.lower().endswith(".csv")
        )
    assert bool(candidates)  # crispron failed?
    if len(candidates) == 1:
        return candidates[0]
    # prefer typical names first
    preferred = [
        p
        for p in candidates
        if os.path.basename(p).lower() in {"crispron.csv", "result.csv", "results.csv"}
    ]
    return preferred[0] if preferred else candidates[0]


def _load_crispron_scores(csv_path: str, expected_30mers: List[str]) -> List[float]:
    """Load CRISPRon scores from a CSV file and align them to expected 30mer guides.

    This function parses the CRISPRon results, matches each scored guide
    to the submitted 30mer sequence, and returns scores ordered by guide index.

    Args:
        csv_path: Path to the CSV file containing CRISPRon output.
        expected_30mers: List of 30mer guide sequences in the order they were 
            submitted.

    Returns:
        A list of CRISPRon scores corresponding to each expected 30mer in input 
            order.

    Raises:
        RuntimeError: If any expected guide is missing a CRISPRon score in the 
            CSV file.
    """
    scores: List[float] = [float("nan")] * len(expected_30mers)
    expected_map = {f"guide_{i}": seq.upper() for i, seq in enumerate(expected_30mers)}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)  # read crispron results
        for row in reader:
            gid = row.get("ID", "")
            seq30 = row.get("30mer", "").upper()
            score_raw = row.get("CRISPRon", "")
            guide_id = gid.split("_")[1]
            if f"guide_{guide_id}" not in expected_map:
                continue
            # keep only the row that matches the submitted 30mer
            if seq30 != expected_map[f"guide_{guide_id}"]:
                continue
            scores[int(guide_id)] = float(score_raw)
    missing = [i for i, s in enumerate(scores) if s != s]  # NaN check
    if missing:
        raise RuntimeError(
            f"Missing CRISPRon scores for {len(missing)} guides. "
            f"First missing indices: {missing[:10]}"
        )
    return scores


def _generate_crispron_tmp_data(tmpdir: str) -> Tuple[str, str]:
    """Create CRISPRon input and output paths inside a temporary directory.

    This function prepares the file system locations that CRISPRon expects
    for the guide FASTA file and the result directory.

    Args:
        tmpdir: Path to the temporary base directory where CRISPRon data will 
            be created.

    Returns:
        A tuple containing the path to the guides FASTA file and the CRISPRon
        results directory within the temporary directory.
    """
    # generate guides fasta required and out folder by crispron script
    return os.path.join(tmpdir, "guides.fa"), os.path.join(tmpdir, "crispron_results")


def _write_guides_fasta(guides: List[str], crispron_fasta: str) -> None:
    """Write guide sequences to a FASTA file in the format expected by CRISPRon.

    This function assigns a unique identifier to each guide and saves them
    in uppercase FASTA format for downstream CRISPRon processing.

    Args:
        guides: List of guide sequences to write to the FASTA file.
        crispron_fasta: Path to the FASTA file where guides will be written.
    """
    with open(crispron_fasta, mode="w") as fout:
        for i, seq in enumerate(guides):
            fout.write(f">guide_{i}\n{seq.upper()}\n")
    assert os.stat(crispron_fasta).st_size > 0


def compute_crispron_score(guides: List[str], conda: str, env_name: str) -> List[float]:
    """Compute CRISPRon on-target activity scores for a list of guide sequences.

    This function runs the external CRISPRon tool in a specified conda
    environment and returns one score per input guide in the original order.

    Args:
        guides: List of guide sequences to be scored by CRISPRon.
        conda: Path to the conda executable used to invoke the CRISPRon environment.
        env_name: Name of the conda environment in which CRISPRon is installed.

    Returns:
        A list of CRISPRon scores corresponding to each guide in the input list.

    Raises:
        subprocess.CalledProcessError: If the CRISPRon subprocess execution fails.
        RuntimeError: If CRISPRon does not produce scores for all provided guides.
    """
    assert bool(guides)  # otherwise we shouldn't be here
    # get path to crispron script
    crispron_root = os.path.abspath(os.path.dirname(__file__))
    crispron_bin = os.path.join(crispron_root, "bin", "CRISPRon.sh")
    # score guides with crispron
    with tempfile.TemporaryDirectory(prefix="crispron_", dir=crispron_root) as tmpdir:           
        # generate guides fasta and output folder required by crispron's script
        crispron_fasta, crispron_tmpdir = _generate_crispron_tmp_data(
            os.path.abspath(tmpdir)
        )
        # write guides to crispron fasta
        _write_guides_fasta(guides, crispron_fasta)
        subprocess.run(
            [
                conda,
                "run",
                "-n",
                env_name,
                "bash",
                crispron_bin,
                crispron_fasta,
                crispron_tmpdir,
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
            cwd=crispron_root,
        )
        csv_path = _find_output_csv(crispron_tmpdir)
        return _load_crispron_scores(csv_path, guides)
