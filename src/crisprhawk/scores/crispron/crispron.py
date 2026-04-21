"""CRISPRon on-target scoring wrapper for CRISPR-HAWK."""

from typing import List
import csv
import tempfile
import subprocess
import os


def _find_output_csv(outdir: str) -> str:
    candidates = []
    for root, _, files in os.walk(outdir):
        for fname in files:
            if fname.lower().endswith(".csv"):
                candidates.append(os.path.join(root, fname))
    if not candidates:
        raise FileNotFoundError(f"No CSV output found in CRISPRon output folder: {outdir}")
    if len(candidates) == 1:
        return candidates[0]

    # prefer typical names first
    preferred = [p for p in candidates if os.path.basename(p).lower() in {"crispron.csv", "result.csv", "results.csv"}]
    return preferred[0] if preferred else candidates[0]


def _load_crispron_scores(csv_path: str, expected_30mers: List[str]) -> List[float]:
    scores: List[float] = [float("nan")] * len(expected_30mers)
    expected_map = {f"guide_{i}": seq.upper() for i, seq in enumerate(expected_30mers)}

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row.get("ID", "")
            seq30 = row.get("30mer", "").upper()
            score_raw = row.get("CRISPRon", "")
            id = gid.split("_")[1]

            if f"guide_{id}" not in expected_map:
                continue

            # keep only the row that matches the submitted 30mer
            if seq30 != expected_map[f"guide_{id}"]:
                continue
            
            scores[int(id)] = float(score_raw)
    missing = [i for i, s in enumerate(scores) if s != s]  # NaN check
    if missing:
        raise RuntimeError(
            f"Missing CRISPRon scores for {len(missing)} guides. "
            f"First missing indices: {missing[:10]}"
        )
    return scores


def compute_crispron_score(
    guides: List[str],
    env_name: str,
    tmp_parent: str,
) -> List[float]:
    """Run CRISPRon on a batch of 30mers and return scores in input order."""
    if not guides:
        return []

    for g in guides:
        if len(g) != 30:
            raise ValueError(
                f"CRISPRon expects 30mer input. Got length {len(g)} for sequence: {g}"
            )

    crispron_root = os.path.abspath(os.path.dirname(__file__))
    crispron_bin = os.path.join(crispron_root, "bin", "CRISPRon.sh")

    if not os.path.isfile(crispron_bin):
        raise FileNotFoundError(f"Cannot find CRISPRon executable: {crispron_bin}")

    with tempfile.TemporaryDirectory(prefix="crispron_", dir=crispron_root) as tmpdir:
        tmpdir = os.path.abspath(tmpdir)
        fasta_path = os.path.abspath(os.path.join(tmpdir, "guides.fa"))
        outdir = os.path.abspath(os.path.join(tmpdir, "crispron_results"))

        with open(fasta_path, "w") as f:
            for i, seq in enumerate(guides):
                f.write(f">guide_{i}\n{seq.upper()}\n")

        subprocess.run(
            ["conda", "run", "-n", env_name, "bash", crispron_bin, fasta_path, outdir],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
            cwd=crispron_root,
        )

        csv_path = _find_output_csv(outdir)
        return _load_crispron_scores(csv_path, guides)