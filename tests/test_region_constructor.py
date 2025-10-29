from crisprhawk import region_constructor

import pytest


class DummyFasta:
    def __init__(self, contig):
        self.contig = contig


class DummyBed:
    def __init__(self, regions):
        self.regions = regions

    def __len__(self):
        return len(self.regions)

    def extract_regions(self, fastas):
        return self.regions


def test_read_fasta(monkeypatch):
    # Patch Fasta to DummyFasta
    monkeypatch.setattr(
        region_constructor,
        "Fasta",
        lambda path, verbosity, debug, faidx=None: DummyFasta(path),
    )
    fastafiles = ["chr1.fa", "chr2.fa"]
    result = region_constructor.read_fasta(fastafiles, 0, False)
    assert isinstance(result, dict)
    assert set(result.keys()) == set(fastafiles)
    for v in result.values():
        assert isinstance(v, DummyFasta)


def test_read_bed(monkeypatch):
    # Patch Bed to DummyBed
    monkeypatch.setattr(
        region_constructor,
        "Bed",
        lambda path, padding, debug: DummyBed(["region1", "region2"]),
    )
    bed = region_constructor.read_bed("regions.bed", 0, False)
    assert isinstance(bed, DummyBed)
    assert len(bed) == 2


def test_extract_regions():
    bed = DummyBed(["regionA", "regionB"])
    fastas = {"chr1": DummyFasta("chr1")}
    regions = region_constructor.extract_regions(bed, fastas, 0, False)
    assert regions == ["regionA", "regionB"]


def test_construct_regions(monkeypatch):
    monkeypatch.setattr(
        region_constructor,
        "read_fasta",
        lambda fastafiles, verbosity, debug: {"chr1": DummyFasta("chr1")},
    )
    monkeypatch.setattr(
        region_constructor,
        "read_bed",
        lambda bedfile, verbosity, debug: DummyBed(["regionX"]),
    )
    monkeypatch.setattr(
        region_constructor,
        "extract_regions",
        lambda bed, fastas, verbosity, debug: ["regionX"],
    )
    regions = region_constructor.construct_regions(
        ["chr1.fa"], "regions.bed", 0, False
    )
    assert regions == ["regionX"]
