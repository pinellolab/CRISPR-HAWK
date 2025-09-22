from crisprhawk.graphical_reports import (
    compute_guide_id,
    assign_guide_type,
    _count_guide_type,
    compute_guide_coord,
    _add_guide_ids,
    _compute_group_delta,
    calculate_deltas,
    _build_alternatives_list,
    _process_single_guide,
    _build_guide_data,
)
from crisprhawk.region import Region

import pandas as pd
import numpy as np

import pytest


def test_compute_guide_id():
    result = compute_guide_id("chr1", 100, 120, "+", "ACGT", "NGG")
    assert result == "chr1_100_120_+_ACGT_NGG"


def test_assign_guide_type_ref():
    assert assign_guide_type("ref", "ACGT", "NGG", False) == 0


def test_assign_guide_type_spacer_pam_alt():
    # Both sgrna and pam are lowercase
    assert assign_guide_type("alt", "acgt", "ngg", False) == 1


def test_assign_guide_type_spacer_alt():
    # sgrna is lowercase, pam is uppercase
    assert assign_guide_type("alt", "acgt", "NGG", False) == 2


def test_assign_guide_type_pam_alt():
    # sgrna is uppercase, pam is lowercase
    assert assign_guide_type("alt", "ACGT", "ngg", False) == 3


def test_count_guide_type():
    guide_types = [0, 1, 2, 3, 0, 1]
    result = _count_guide_type(guide_types)
    assert result["Reference Guides"] == 2
    assert result["Spacer+PAM Alternative Guides"] == 2
    assert result["Spacer Alternative Guides"] == 1
    assert result["PAM Alternative Guides"] == 1


def test_compute_guide_coord():
    result = compute_guide_coord("chr2", 200, 220, "-")
    assert result == "chr2_200_220_-"


def test_add_guide_ids():
    df = pd.DataFrame(
        [
            ["chr1", 100, 120, "foo", "bar", "baz", "+"],
            ["chr2", 200, 220, "foo", "bar", "baz", "-"],
        ]
    )
    df2 = _add_guide_ids(df, False)
    assert "guide_id" in df2.columns
    assert df2["guide_id"].iloc[0] == "chr1_100_120_+"
    assert df2["guide_id"].iloc[1] == "chr2_200_220_-"


def test_compute_group_delta():
    df = pd.DataFrame({"origin": ["ref", "alt", "alt"], "score": [0.8, 0.5, 0.6]})
    result = _compute_group_delta(df, "score")
    assert np.isclose(result["delta"].iloc[0], 0.0)
    assert np.isclose(result["delta"].iloc[1], -0.3)
    assert np.isclose(result["delta"].iloc[2], -0.2)


def test_calculate_deltas():
    df = pd.DataFrame(
        {
            "guide_id": ["g1", "g1", "g2"],
            "origin": ["ref", "alt", "ref"],
            "score": [1.0, 0.7, 0.5],
        }
    )
    result = calculate_deltas(df, "score")
    assert "delta" in result.columns
    assert np.isclose(result[result["guide_id"] == "g1"]["delta"].iloc[0], 0.0)
    assert np.isclose(result[result["guide_id"] == "g1"]["delta"].iloc[1], -0.3)
    assert np.isclose(result[result["guide_id"] == "g2"]["delta"].iloc[0], 0.0)


def test_build_alternatives_list():
    df = pd.DataFrame(
        {
            "sgRNA_sequence": ["acgt", "tgca"],
            "pam": ["ngg", "ncc"],
            "score_col": [0.5, 0.4],
            "delta": [-0.3, -0.4],
            "samples": [10, 20],
            "variant_id": ["v1", "v2"],
        }
    )
    result = _build_alternatives_list(df, "score_col")
    assert result[0]["alt_sgRNA"] == "acgt"
    assert result[1]["pam"] == "ncc"
    assert result[0]["alt_score"] == 0.5


def test_process_single_guide():
    df = pd.DataFrame(
        {
            "origin": ["ref", "alt", "alt"],
            "sgRNA_sequence": ["ACGT", "acgt", "tgca"],
            "pam": ["NGG", "ngg", "ncc"],
            "score_col": [0.8, 0.5, 0.6],
            "delta": [0.0, -0.3, -0.2],
            "samples": [10, 20, 30],
            "variant_id": ["v0", "v1", "v2"],
        }
    )
    result = _process_single_guide(df, "score_col")
    assert result["ref_sgRNA"] == "ACGT"
    assert result["pam"] == "NGG"
    assert result["ref_score"] == 0.8
    assert len(result["alts"]) == 2


def test_build_guide_data():
    df = pd.DataFrame(
        {
            "guide_id": ["g1", "g1", "g2"],
            "origin": ["ref", "alt", "ref"],
            "sgRNA_sequence": ["ACGT", "acgt", "TGCA"],
            "pam": ["NGG", "ngg", "NCC"],
            "score_col": [0.8, 0.5, 0.6],
            "delta": [0.0, -0.3, 0.0],
            "samples": [10, 20, 30],
            "variant_id": ["v0", "v1", "v2"],
        }
    )
    result = _build_guide_data(df, "score_col")
    assert "g1" in result
    assert "g2" in result
    assert result["g1"]["ref_sgRNA"] == "ACGT"
    assert result["g2"]["ref_sgRNA"] == "TGCA"
