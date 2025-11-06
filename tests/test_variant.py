import pytest
from crisprhawk import variant


def test_assign_vtype():
    assert variant._assign_vtype("A", "T") == "snp"
    assert variant._assign_vtype("A", "AT") == "indel"
    assert variant._assign_vtype("AT", "A") == "indel"


def test_compute_id():
    assert variant._compute_id("chr1", 123, "A", "T") == "chr1-123-A/T"
    assert variant._compute_id("2", 456, "G", "C") == "2-456-G/C"


def test_adjust_multiallelic_snp():
    ref, alt, pos = variant.adjust_multiallelic("A", "T", 100)
    assert (ref, alt, pos) == ("A", "T", 100)


def test_adjust_multiallelic_deletion():
    ref, alt, pos = variant.adjust_multiallelic("AAT", "AA", 200)
    assert (ref, alt, pos) == ("AT", "A", 201)


def test_adjust_multiallelic_insertion():
    ref, alt, pos = variant.adjust_multiallelic("AA", "AAT", 300)
    assert (ref, alt, pos) == ("A", "AT", 301)


def test_variantrecord_repr_and_str():
    v = variant.VariantRecord(debug=False)
    v._chrom = "chr1"
    v._position = 123
    v._ref = "A"
    v._alt = ["T"]
    assert "chr1" in repr(v)
    assert "A" in str(v)
    assert "T" in str(v)


def test_variantrecord_eq_and_hash():
    v1 = variant.VariantRecord(debug=False)
    v1._chrom = "chr1"
    v1._position = 123
    v1._ref = "A"
    v1._alt = ["T"]
    v2 = variant.VariantRecord(debug=False)
    v2._chrom = "chr1"
    v2._position = 123
    v2._ref = "A"
    v2._alt = ["T"]
    v2.contig  # property access
    v2.position
    v2.ref
    v2.alt
    assert v1 == v2
    assert hash(v1) == hash(v2)


def test_variantrecord_split_and_get_altalleles():
    v = variant.VariantRecord(debug=False)
    v._chrom = "chr1"
    v._position = 123
    v._ref = "A"
    v._alt = ["T", "AT"]
    v._allelesnum = 2
    v._vtype = ["snp", "indel"]
    v._filter = "PASS"
    v._vid = ["id1", "id2"]
    v._afs = [0.1, 0.2]
    v._samples = [({"sample1"}, set()), (set(), {"sample2"})]
    snps = v.split("snp")
    indels = v.split("indel")
    assert len(snps) == 1
    assert snps[0].alt == ["T"]
    assert len(indels) == 1
    assert indels[0].alt == ["AT"]
    assert v.get_altalleles("snp") == ["T"]
    assert v.get_altalleles("indel") == ["AT"]


def test_variantrecord_pytest_initialize():
    v = variant.VariantRecord(debug=False)
    v.pytest_initialize(100, "A", "T", "snp", "id1", [0.5])
    assert v._chrom == "chrx"
    assert v._position == 100
    assert v._ref == "A"
    assert v._alt == ["T"]
    assert v._vtype == ["snp"]
    assert v._vid == ["id1"]
    assert v._afs == [0.5]


def test_retrieve_alt_alleles():
    v = variant.VariantRecord(debug=False)
    alleles = v._retrieve_alt_alleles("A,T,G")
    assert alleles == ["A", "T", "G"]
    assert v._retrieve_alt_alleles("") == []


def test_find_tbi(tmp_path):
    vcf_path = tmp_path / "test.vcf"
    tbi_path = tmp_path / "test.vcf.tbi"
    vcf_path.write_text("dummy")
    assert not variant.find_tbi(str(vcf_path))
    tbi_path.write_text("dummy")
    assert variant.find_tbi(str(vcf_path))
