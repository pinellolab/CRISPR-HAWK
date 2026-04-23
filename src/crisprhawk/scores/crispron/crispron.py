""" """

from typing import List, Tuple

import csv
import os
import subprocess
import tempfile


def _find_output_csv(crispron_tmpdir: str) -> str:
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
    # generate guides fasta required and out folder by crispron script
    return os.path.join(tmpdir, "guides.fa"), os.path.join(tmpdir, "crispron_results")


def _write_guides_fasta(guides: List[str], crispron_fasta: str) -> None:
    with open(crispron_fasta, mode="w") as fout:
        for i, seq in enumerate(guides):
            fout.write(f">guide_{i}\n{seq.upper()}\n")
    assert os.stat(crispron_fasta).st_size > 0


def compute_crispron_score(guides: List[str], conda: str, env_name: str) -> List[float]:
    """Run CRISPRon on a batch of 30mers and return scores in input order."""
    assert bool(guides)  # otherwise we shouldn't be here
    # get path to crispron script
    crispron_root = os.path.abspath(os.path.dirname(__file__))
    crispron_bin = os.path.join(crispron_root, "bin", "CRISPRon.sh")
    
    # score guides with crispron
    # with tempfile.TemporaryDirectory(prefix="crispron_", dir=crispron_root) as tmpdir:
    tmpdir = tempfile.mkdtemp(prefix="crispron_", dir=crispron_root)
        
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
