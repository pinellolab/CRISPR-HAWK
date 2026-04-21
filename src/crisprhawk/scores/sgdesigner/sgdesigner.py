"""sgDesigner scoring wrapper for CRISPR-HAWK."""

from typing import List
import csv
import os
import subprocess
import tempfile


def _find_output_txt(outdir: str) -> str:
    candidates = []
    for root, _, files in os.walk(outdir):
        for fname in files:
            if fname == "sgDesigner_V2.0_prediction_result.txt":
                candidates.append(os.path.join(root, fname))

    if not candidates:
        raise FileNotFoundError(
            f"No sgDesigner output file found in output folder: {outdir}"
        )

    return candidates[0]


def _load_sgdesigner_scores(txt_path: str, expected_26mers: List[str]) -> List[float]:
    """
    Parse sgDesigner output and return one score per input sequence in input order.

    Expected output format:
        seqId    Score    Sequence    Orientation    Position

    We keep only rows where:
    - seqId corresponds to one of our submitted FASTA headers (guide_i)
    - Sequence matches the first 20 nt of the submitted 26mer
      because sgDesigner reports only the spacer sequence in the output file
    """
    scores: List[float] = [float("nan")] * len(expected_26mers)
    expected_map = {f"guide_{i}": seq[:20].upper() for i, seq in enumerate(expected_26mers)}
    with open(txt_path, newline="\n") as f:
        all_lines = f.readlines()
        for row in all_lines[1:]:
            gid, score_raw, spacer, _, _ = row.split("\t")
            # gid = row.get("seqId", "")
            # spacer = row.get("Sequence", "").upper()
            # score_raw = row.get("Score", "")
            id = gid.split("_")[1]

            if f"guide_{id}" not in expected_map:
                continue

            # keep only the row matching the submitted 20-nt spacer
            if spacer.upper() != expected_map[f"guide_{id}"]:
                continue

            scores[int(id)] = float(score_raw)
    missing = [i for i, s in enumerate(scores) if s != s]
    if missing:
        raise RuntimeError(
            f"Missing sgDesigner scores for {len(missing)} guides. "
            f"First missing indices: {missing[:10]}"
        )

    return scores


def compute_sgdesigner_score(
    guides: List[str],
    env_name: str,
    tmp_parent: str,
) -> List[float]:
    """
    Run sgDesigner on a batch of 26mers and return scores in input order.

    Input expected by sgDesigner wrapper:
        20 nt spacer + 3 nt PAM + 3 nt downstream 
    """
    if not guides:
        return []

    for g in guides:
        if len(g) != 26:
            raise ValueError(
                f"sgDesigner expects 26mer input. Got length {len(g)} for sequence: {g}"
            )

    sgdesigner_root = os.path.abspath(os.path.dirname(__file__))
    sgdesigner_pl = os.path.join(sgdesigner_root, "sgDesigner.pl")

    if not os.path.isfile(sgdesigner_pl):
        raise FileNotFoundError(f"Cannot find sgDesigner script: {sgdesigner_pl}")


    with tempfile.TemporaryDirectory(prefix="sgdesigner_", dir=sgdesigner_root) as tmpdir:
        tmpdir = os.path.abspath(tmpdir)
        fasta_path = os.path.join(tmpdir, "guides.fa")
        result_dir = os.path.join(tmpdir, "result")
        temp_dir = os.path.join(tmpdir, "temp")

        os.makedirs(result_dir, exist_ok=True)
        os.makedirs(temp_dir, exist_ok=True)

        with open(fasta_path, "w") as f:
            for i, seq in enumerate(guides):
                f.write(f">guide_{i}\n{seq.upper()}\n")

        env = os.environ.copy()
        env["SGDESIGNER_RESULT_DIR"] = result_dir
        env["SGDESIGNER_TEMP_DIR"] = temp_dir

        subprocess.run(
            ["conda", "run", "-n", env_name, "perl", sgdesigner_pl, "-f", fasta_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
            cwd=sgdesigner_root,
            env=env
        )
        txt_path = _find_output_txt(result_dir)
        return _load_sgdesigner_scores(txt_path, guides)