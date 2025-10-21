from crisprhawk.region import Region, RegionList, RegionListIterator
from crisprhawk.coordinate import Coordinate
from crisprhawk.sequence import Sequence

import pytest


def make_region(contig="chr1", start=1, stop=10, seq="ACTGACTGAC"):
    return Region(Sequence(seq, False), Coordinate(contig, start, stop, 0))


def test_region_len():
    region = make_region(seq="ACTG")
    assert len(region) == 4


def test_region_eq_and_hash():
    r1 = make_region()
    r2 = make_region()
    r3 = make_region(seq="TTTT")
    assert r1 == r2
    assert r1 != r3
    assert hash(r1) == hash(r2)
    assert hash(r1) != hash(r3)


def test_region_str_and_repr():
    region = make_region()
    assert str(region).startswith(">")
    assert "chr1" in str(region)
    assert "ACTGACTGAC" in str(region)
    assert "<Region object;" in repr(region)


def test_region_getitem():
    region = make_region(seq="ACTGACTGAC")
    assert region[0] == "A"
    assert region[1:4] == list("CTG")


def test_region_contain_overlap():
    r1 = make_region(start=1, stop=10)
    r2 = make_region(start=2, stop=5)
    r3 = make_region(start=6, stop=15)
    assert r1.contains(r2)
    assert not r2.contains(r1)
    assert r1.overlap(r2)
    assert r1.overlap(r3)
    assert not r2.overlap(r3)


def test_region_properties():
    region = make_region(contig="chr2", start=5, stop=15, seq="TTTTTTTTTT")
    assert region.contig == "chr2"
    assert region.start == 5
    assert region.stop == 15
    assert isinstance(region.sequence, Sequence)
    assert isinstance(region.coordinates, Coordinate)


def test_regionlist_len_and_iter():
    regions = [make_region(seq="A"), make_region(seq="T")]
    region_list = RegionList(regions)
    assert len(region_list) == 2
    assert list(region_list) == regions


def test_regionlist_repr_and_str():
    regions = [make_region(seq="A"), make_region(seq="T")]
    region_list = RegionList(regions)
    assert "<RegionList object;" in repr(region_list)
    assert "chr1" in str(region_list)


def test_regionlist_getitem():
    regions = [make_region(seq="A"), make_region(seq="T")]
    region_list = RegionList(regions)
    assert region_list[0] == regions[0]
    assert region_list[1] == regions[1]
    with pytest.raises(IndexError):
        _ = region_list[2]


def test_regionlist_extend_and_append():
    r1 = make_region(seq="A")
    r2 = make_region(seq="T")
    region_list1 = RegionList([r1])
    region_list2 = RegionList([r2])
    region_list1.extend(region_list2)
    assert len(region_list1) == 2
    r3 = make_region(seq="G")
    region_list1.append(r3)
    assert len(region_list1) == 3
    with pytest.raises(TypeError):
        region_list1.append("not_a_region") # type: ignore
    with pytest.raises(TypeError):
        region_list1.extend(["not_a_regionlist"]) # type: ignore


def test_regionlistiterator_next():
    regions = [make_region(seq="A"), make_region(seq="T")]
    region_list = RegionList(regions)
    iterator = RegionListIterator(region_list)
    assert next(iterator) == regions[0]
    assert next(iterator) == regions[1]
    with pytest.raises(StopIteration):
        next(iterator)
