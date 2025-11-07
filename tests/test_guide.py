from crisprhawk.guide import Guide, GuideIterator
from crisprhawk.crisprhawk_error import CrisprHawkGuideError

import numpy as np

import pytest


class DummyException(Exception):
    pass


def test_guide_initialization():
    guide = Guide(
        position_start=100,
        position_stop=123,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={"var1": 0.5},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    assert guide.start == 100
    assert guide.stop == 123
    assert guide.sequence.startswith("N" * 10)
    assert guide.samples == "sample1"
    assert guide.variants == "var1"
    assert guide.afs == {"var1": 0.5}
    assert guide.hapid == "hap1"
    assert isinstance(guide.guide_id, str)
    assert guide.azimuth_score == "NA"
    assert guide.rs3_score == "NA"
    assert guide.deepcpf1_score == "NA"
    assert guide.cfdon_score == "NA"
    assert guide.elevationon_score == "NA"
    assert guide.ooframe_score == "NA"
    assert guide.gc == "NA"
    assert guide.funcann == []
    assert guide.geneann == []
    assert guide.offtargets == "NA"


def test_guide_repr_and_len():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    rep = repr(guide)
    assert "start=1" in rep
    assert "stop=24" in rep
    assert len(guide) == len(guide.sequence)


def test_guide_getitem_and_iter():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    # __getitem__
    assert guide[0] == "N"
    assert guide[:4] == "NNNN"
    # __iter__ and GuideIterator
    seq = "".join([nt for nt in guide])
    assert seq == guide.sequence


def test_guide_reverse_complement():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    orig_seq = guide.sequence
    guide.reverse_complement()
    assert guide.sequence != orig_seq
    assert isinstance(guide.sequence, str)


def test_guide_setters_and_properties():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    guide.variants = "var2"
    assert guide.variants == "var2"
    guide.afs_str = ["0.1", "0.2"]
    assert guide.afs_str == "0.1,0.2"
    guide.azimuth_score = 0.5
    assert guide.azimuth_score == "0.5"
    guide.rs3_score = 0.7
    assert guide.rs3_score == "0.7"
    guide.deepcpf1_score = 0.8
    assert guide.deepcpf1_score == "0.8"
    guide.cfdon_score = 0.9
    assert guide.cfdon_score == "0.9"
    guide.elevationon_score = 1.0
    assert guide.elevationon_score == "1.0"
    guide.ooframe_score = 2
    assert guide.ooframe_score == "2"
    guide.gc = 0.45
    assert guide.gc == "0.45"
    guide.funcann = "ann1"
    assert "ann1" in guide.funcann
    guide.geneann = "gene1"
    assert "gene1" in guide.geneann
    guide.offtargets = 3
    assert guide.offtargets == "3"


def test_guide_geneann_type_error():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    with pytest.raises(CrisprHawkGuideError):
        guide.geneann = 123  # not a string


def test_guide_offtargets_type_error():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    with pytest.raises(CrisprHawkGuideError):
        guide.offtargets = "not_an_int"


def test_guide_cfd_type_error():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    with pytest.raises(CrisprHawkGuideError):
        guide.cfd = "not_a_float"


def test_guideiterator_next_and_stopiteration():
    guide = Guide(
        position_start=1,
        position_stop=24,
        sequence="N" * 10 + "AGCTTAGCTAGCTAGCTAGCTAG" + "N" * 10,
        guidelen=20,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        posmap={i:i for i in range(10 + 20 + 3 + 10)},
        debug=False,
        right=True,
        hapid="hap1",
    )
    iterator = GuideIterator(guide)
    collected = []
    try:
        while True:
            collected.append(next(iterator))
    except StopIteration:
        pass
    assert "".join(collected) == guide.sequence
