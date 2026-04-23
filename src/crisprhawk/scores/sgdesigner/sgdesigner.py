""" """

from ...utils import create_folder

from typing import List, Tuple, Dict

import os
import subprocess
import tempfile


def _find_output_txt(sgdesigner_outdir: str) -> str:
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
    scores: List[float] = [float("nan")] * len(expected_26mers)
    expected_map = {f"guide_{i}": seq[:20].upper() for i, seq in enumerate(expected_26mers)}
    with open(txt_path, mode="r") as fin:
        fin.readline()  # skip header
        for line in fin:
            gid, score_raw, spacer, _, _ = line.split("\t")
            guide_id = gid.split("_")[1]
            if f"guide_{id}" not in expected_map:
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
    # generate guides fasta required and out folder by crispron script
    return os.path.join(tmpdir, "guides.fa"), os.path.join(tmpdir, "sgdesigner_results"), os.path.join(tmpdir, "temp")

def _write_guides_fasta(guides: List[str], sgdesigner_fasta: str) -> None:
    with open(sgdesigner_fasta, mode="w") as fout:
        for i, seq in enumerate(guides):
            fout.write(f">guide_{i}\n{seq.upper()}\n")
    assert os.stat(sgdesigner_fasta).st_size > 0

def _init_environ(sgdesigner_results: str, sgdesigner_tmp: str) -> Dict[str, str]:
    env = os.environ.copy()
    env["SGDESIGNER_RESULT_DIR"] = sgdesigner_results
    env["SGDESIGNER_TEMP_DIR"] = sgdesigner_tmp
    return env 


def compute_sgdesigner_score(guides: List[str], env_name: str) -> List[float]:
    assert bool(guides)  # otherwise we shouldn't be here
    # get path to sgdesigner script
    sgdesigner_root = os.path.abspath(os.path.dirname(__file__))
    sgdesigner_pl = os.path.join(sgdesigner_root, "sgDesigner.pl")
    # score guides with sgdesigner
    with tempfile.TemporaryDirectory(prefix="sgdesigner_", dir=sgdesigner_root) as tmpdir:
        # generate guides fasta, output folder, and tmp folder required 
        # by sgdesigner's script
        sgdesigner_fasta, sgdesigner_results, sgdesigner_tmpdir = _generate_sgdesigner_tmp_data(tmpdir)
        for d in [sgdesigner_results, sgdesigner_tmpdir]:
            create_folder(d, exist_ok=True)  # create sgdesigner folders
        # write guides to sgdesigner fasta
        _write_guides_fasta(guides, sgdesigner_fasta)

        subprocess.run(
            ["conda", "run", "-n", env_name, "perl", sgdesigner_pl, "-f", sgdesigner_fasta],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
            cwd=sgdesigner_root,
            env=_init_environ(sgdesigner_results, sgdesigner_tmpdir)
        )
        txt_path = _find_output_txt(sgdesigner_results)
        return _load_sgdesigner_scores(txt_path, guides)