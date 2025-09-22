from crisprhawk.reports import (
    compute_pam_class,
    compute_guide_origin,
    compute_strand_orientation,
    collapse_samples,
    parse_variant_ids,
    check_variant_ids,
    collapse_haplotype_ids,
    collapse_annotation,
    collapse_offtargets,
    collapse_cfd,
)
from crisprhawk.pam import PAM
from crisprhawk.region import Region
from crisprhawk.guide import Guide

import pandas as pd

import pytest


def test_compute_pam_class():
    pam = PAM("NGG", False, False)
    assert compute_pam_class(pam) == "[ACGT]GG"
    pam2 = PAM("TTTV", True, False)
    assert compute_pam_class(pam2) == "TTT[ACG]"


def test_compute_guide_origin():
    assert compute_guide_origin("REF") == "ref"
    assert compute_guide_origin("ALT") == "alt"


def test_compute_strand_orientation():
    assert compute_strand_orientation(0) == "+"
    assert compute_strand_orientation(1) == "-"


def test_collapse_samples():
    s = pd.Series(["sample1:1|0", "sample2:0|1"])
    result = collapse_samples(s)
    assert "sample1:1|0" in result or "sample2:0|1" in result


def test_parse_variant_ids():
    assert parse_variant_ids("var1,var2") == {"var1", "var2"}
    assert parse_variant_ids("") == set()


def test_check_variant_ids():
    assert (
        check_variant_ids(["var1,var2", "var2,var1"]) == "var1,var2"
        or check_variant_ids(["var1,var2", "var2,var1"]) == "var2,var1"
    )


def test_collapse_haplotype_ids():
    s = pd.Series(["hap1,hap2", "hap2,hap3"])
    result = collapse_haplotype_ids(s)
    assert set(result.split(",")) == {"hap1", "hap2", "hap3"}


def test_collapse_annotation():
    s = pd.Series(["ann1,ann2", "ann2,ann3"])
    result = collapse_annotation(s)
    assert set(result.split(",")) == {"ann1", "ann2", "ann3"}


def test_collapse_offtargets():
    s = pd.Series([5, 5, 5])
    assert collapse_offtargets(s) == 5


def test_collapse_cfd():
    s = pd.Series(["0.5", "0.5"])
    assert collapse_cfd(s) == "0.5"
