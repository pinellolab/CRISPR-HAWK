from crisprhawk.haplotypes import (
    initialize_haplotypes,
    ishomozygous,
    collapse_haplotypes,
    classify_variants,
    generate_haplotype_ids,
)

import pytest


class DummySequence:
    def __init__(self, sequence, debug, allow_lower_case=False):
        self.sequence = sequence


class DummyCoordinate:
    def __init__(self, contig, start, stop, strand):
        self.contig = contig
        self.start = start
        self.stop = stop
        self.startp = start
        self.stopp = stop


class DummyRegion:
    def __init__(self, contig, start, stop, seq):
        self.contig = contig
        self.start = start
        self.stop = stop
        self.coordinates = DummyCoordinate(contig, start, stop, 0)
        self.sequence = DummySequence(seq, False)


class DummyVariantRecord:
    def __init__(self, vtype, samples, position=0, ref="A", alt=None, id=None):
        self.vtype = (vtype,)
        self.samples = samples
        self.position = position
        self.ref = ref
        self.alt = alt or []
        self.id = id or ["1-1-A/T"]


class DummyHaplotype:
    def __init__(self, sequence, coordinates, phased, copy, debug):
        self.sequence = sequence
        self.coordinates = coordinates
        self.phased = phased
        self.copy = copy
        self.debug = debug
        self.samples = ""
        self.variants = ""
        self.id = None
        self.afs = None
        self.posmap = None
        self.posmap_rev = None
        self.variant_alleles = None

    def add_variants_phased(self, variants, sample):
        pass

    def add_variants_unphased(self, variants, sample):
        pass

    def homozygous_samples(self):
        self.samples = "homozygous"

    def set_afs(self, afs):
        self.afs = afs

    def set_posmap(self, posmap, posmap_rev):
        self.posmap = posmap
        self.posmap_rev = posmap_rev

    def set_variant_alleles(self, alleles):
        self.variant_alleles = alleles


@pytest.fixture
def dummy_region():
    return DummyRegion("chr1", 0, 100, "A" * 101)


def test_initialize_haplotypes(dummy_region, monkeypatch):
    # Patch Haplotype to DummyHaplotype for isolated testing
    import crisprhawk.haplotypes as haplomod

    monkeypatch.setattr(haplomod, "Haplotype", DummyHaplotype)
    monkeypatch.setattr(haplomod, "Sequence", DummySequence)
    regions = [dummy_region]
    result = initialize_haplotypes(regions, debug=False)
    assert len(result) == 1
    assert isinstance(list(result.values())[0][0], DummyHaplotype)


def test_ishomozygous():
    h1 = DummyHaplotype(DummySequence("AAAA", False), None, False, 0, False)
    h2 = DummyHaplotype(DummySequence("AAAA", False), None, False, 1, False)
    assert ishomozygous([h1, h2]) is True
    h3 = DummyHaplotype(DummySequence("AAAT", False), None, False, 1, False)
    assert ishomozygous([h1, h3]) is False


def test_collapse_haplotypes(monkeypatch):
    import crisprhawk.haplotypes as haplomod

    monkeypatch.setattr(haplomod, "Haplotype", DummyHaplotype)
    monkeypatch.setattr(haplomod, "Sequence", DummySequence)
    h1 = DummyHaplotype(DummySequence("AAAA", False), None, False, 0, False)
    h2 = DummyHaplotype(DummySequence("AAAA", False), None, False, 1, False)
    h3 = DummyHaplotype(DummySequence("AAAT", False), None, False, 1, False)
    collapsed = collapse_haplotypes([h1, h2, h3], debug=False)
    assert any(h.sequence.sequence == "AAAA" for h in collapsed)
    assert any(h.sequence.sequence == "AAAT" for h in collapsed)
    assert len(collapsed) == 2


def test_classify_variants():
    v1 = DummyVariantRecord("snp", [[["sample1"]]])
    v2 = DummyVariantRecord("indel", [[["sample2"]]])
    snvs, indels = classify_variants([v1, v2])
    assert v1 in snvs
    assert v2 in indels


def test_generate_haplotype_ids():
    h1 = DummyHaplotype(DummySequence("AAAA", False), None, False, 0, False)
    h2 = DummyHaplotype(DummySequence("AAAT", False), None, False, 1, False)
    region = DummyRegion("chr1", 0, 100, "A" * 101)
    haplotypes = {region: [h1, h2]}
    result = generate_haplotype_ids(haplotypes)
    ids = [h.id for h in result[region]]
    assert len(ids) == 2
    assert all(isinstance(i, str) and i.startswith("hap_") for i in ids)
