from crisprhawk.haplotype import (
    Haplotype,
    HaplotypeIndel,
    _sort_variants,
    _compute_chains,
    _encode_iupac,
)
from crisprhawk.sequence import Sequence
from crisprhawk.coordinate import Coordinate
from crisprhawk.variant import VariantRecord, VTYPES
from crisprhawk.crisprhawk_error import CrisprHawkHaplotypeError, CrisprHawkIupacTableError

import pytest


def make_variant(position, ref, alt, vtype, vid="v1", afs=[0.5]):
    v = VariantRecord(False)
    v.pytest_initialize(position, ref, alt, vtype, vid, afs)
    return v

def test_haplotype_init_and_str():
    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hap = Haplotype(seq, coord, phased=True, chromcopy=0, debug=False)
    assert hap._size == 4
    assert hap._phased is True
    assert str(hap) == "REF: ACGT"

def test_haplotype_setters_and_getters():
    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hap = Haplotype(seq, coord, phased=False, chromcopy=1, debug=False)
    hap.samples = "sample1"
    assert hap.samples == "sample1"
    hap.variants = "v1"
    assert hap.variants == "v1"
    hap.set_afs({"v1": 0.5})
    assert hap.afs == {"v1": 0.5}
    hap.id = "hap1"
    assert hap.id == "hap1"
    hap.set_posmap({0: 0, 1: 1}, {0: 0, 1: 1})
    assert hap.posmap[0] == 0
    hap.set_variant_alleles({0: ("A", "G")})
    assert hap.variant_alleles[0] == ("A", "G")

def test_haplotype_add_variants_phased_and_unphased():
    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hap = Haplotype(seq, coord, phased=True, chromcopy=0, debug=False)
    var1 = make_variant(0, "A", "G", VTYPES[0], vid="v1")
    var2 = make_variant(1, "C", "T", VTYPES[0], vid="v2")
    hap.add_variants_phased([var1, var2], "sample1")
    assert "sample1" in hap.samples
    assert "v1" in hap.variants and "v2" in hap.variants
    assert hap.afs["v1"] == 0.5

    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hap2 = Haplotype(seq, coord, phased=False, chromcopy=1, debug=False)
    hap2.add_variants_unphased([var1, var2], "sample2")
    assert hap2.samples == "sample2"
    assert "v1" in hap2.variants and "v2" in hap2.variants

def test_haplotype_homozygous_samples():
    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hap = Haplotype(seq, coord, phased=True, chromcopy=0, debug=False)
    hap.samples = "sample1:1|0"
    hap.homozygous_samples()
    assert hap.samples == "sample1:1|1"

    hap2 = Haplotype(seq, coord, phased=True, chromcopy=0, debug=False)
    hap2.samples = "REF"
    with pytest.raises(SystemExit):
        hap2.homozygous_samples()

def test_sort_variants_and_compute_chains():
    var1 = make_variant(0, "A", "G", VTYPES[0], vid="v1")
    var2 = make_variant(1, "C", "T", VTYPES[0], vid="v2")
    var3 = make_variant(2, "G", "GA", VTYPES[1], vid="v3")
    sorted_vars = _sort_variants([var3, var1, var2])
    assert sorted_vars[0].id[0] == "v1"
    assert sorted_vars[1].id[0] == "v2"
    assert sorted_vars[2].id[0] == "v3"
    chains = _compute_chains([var1, var3])
    assert chains == [0, 1]

def test_encode_iupac_success_and_failure():
    # A/G = R, C/T = Y, etc.
    assert _encode_iupac("A", "G", 0, False) == "R"
    assert _encode_iupac("C", "T", 1, False) == "Y"
    with pytest.raises(SystemExit):
        _encode_iupac("A", "B", 2, False)

def test_haplotypeindel_properties_and_setters():
    seq = Sequence("ACGT", debug=False)
    coord = Coordinate("chr1", 0, 4, 0)
    hapindel = HaplotypeIndel(seq, coord, phased=True, chromcopy=0, debug=False)
    hapindel.offset = 5
    assert hapindel.offset == 5
    hapindel.indel_position = 2
    assert hapindel.indel_position == 2
    with pytest.raises(Exception):
        hapindel.offset = "bad"
    with pytest.raises(Exception):
        hapindel.indel_position = "bad"