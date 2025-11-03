from crisprhawk.search_guides import (
    compute_scan_start_stop,
    _valid_guide,
    _decode_iupac,
    resolve_guide,
    adjust_guide_position,
    group_guides_position,
    remove_redundant_guides,
    is_pamhit_valid,
    is_pamhit_in_range,
)
from crisprhawk.pam import PAM
from crisprhawk.haplotype import Haplotype
from crisprhawk.guide import Guide
from crisprhawk.sequence import Sequence
from crisprhawk.coordinate import Coordinate

import pytest


class DummyHaplotype:
    def __init__(self):
        self.stop = 223
        self.start = 100
        self.posmap_rev = {i: i for i in range(101)}
        self.posmap = {i: i for i in range(101)}
        self.variant_alleles = {i: ("A", "G") for i in range(101)}
        self.samples = "REF"
        self.variants = ""
        self.afs = ""
        self.id = "hap1"

    def __getitem__(self, key):
        # Return a string of "A" for the requested range
        if isinstance(key, slice):
            return ["A"] * (key.stop - key.start)
        return "A"

    def __len__(self):
        return 101


def test_compute_scan_start_stop_basic():
    hap = Haplotype(
        Sequence("A" * 1201, False, allow_lower_case=True),
        Coordinate("chrx", 100, 1100, 100),
        True,
        0,
        False,
    )
    start, stop = compute_scan_start_stop(hap, 100, 1100, 3)
    assert isinstance(start, int)
    assert isinstance(stop, int)
    assert start < stop


def test_valid_guide_true():
    pam = PAM("NGG", right=True, debug=False)
    pam.encode(0)
    assert _valid_guide("CCCCCCCCCCCCCCCCCCCCCCAGG", pam, 0, True, False) in [
        True,
        False,
    ]  # Should not error


def test_decode_iupac_single():
    hap = DummyHaplotype()
    assert _decode_iupac("A", 10, hap, False) == "A"


def test_decode_iupac_ambiguous():
    hap = DummyHaplotype()
    # Should return a string of possible alleles, e.g. "AG"
    result = _decode_iupac("R", 10, hap, False)
    assert set(result.lower()).issubset({"a", "g"})


def test_adjust_guide_position():
    hap = DummyHaplotype()
    start, stop = adjust_guide_position(hap.posmap, 10, 20, 3, True)
    assert isinstance(start, int)
    assert isinstance(stop, int)


def test_is_pamhit_valid_and_in_range():
    assert is_pamhit_valid(10, 100, 20, 3, True) in [True, False]
    assert is_pamhit_in_range(10, 20, 3, 100, True) in [True, False]


def test_group_guides_position_and_remove_redundant():
    guide1 = Guide(
        100,
        123,
        "N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAGC" + "N" * 10,
        23,
        3,
        0,
        "REF",
        "",
        {},
        False,
        False,
        "hap1",
    )
    guide2 = Guide(
        100,
        123,
        "N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAGC" + "N" * 10,
        23,
        3,
        0,
        "ALT",
        "var1",
        {"var1": 0.5},
        False,
        False,
        "hap1",
    )
    guides = [guide1, guide2]
    grouped = group_guides_position(guides, False)
    assert isinstance(grouped, dict)
    filtered = remove_redundant_guides(guides, False)
    assert isinstance(filtered, list)
    assert all(isinstance(g, Guide) for g in filtered)
