from crisprhawk.scoring import (
    azimuth_score,
    rs3_score,
    deepcpf1_score,
    group_guides_position,
    cfdon_score,
    elevationon_score,
    outofframe_score,
)
from crisprhawk.guide import Guide
from crisprhawk.utils import prepare_package

import pytest



def make_guide(**kwargs):
    defaults = dict(
        position_start=1,
        position_stop=24,
        sequence="N" * 50 + "AGCTTAGCTAGCTAGCTAGCTAGC" + "N" * 50,
        guidelen=23,
        pamlen=3,
        direction=1,
        samples="sample1",
        variants="var1",
        afs={},
        debug=False,
        right=True,
        hapid="hap1",
    )
    defaults.update(kwargs)
    return Guide(**defaults)


def test_azimuth_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.scoring.azimuth", lambda seqs: [0.5])
    scored = azimuth_score(guides, verbosity=0, debug=False)
    assert scored[0].azimuth_score == "0.5"


def test_rs3_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.scoring.rs3", lambda seqs: [0.7])
    scored = rs3_score(guides, verbosity=0, debug=False)
    assert scored[0].rs3_score == "0.7"


def test_deepcpf1_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.scoring.deepcpf1", lambda seqs: [0.8])
    scored = deepcpf1_score(guides, verbosity=0, debug=False)
    assert scored[0].deepcpf1_score == "0.8"


def test_group_guides_position():
    guides = [make_guide(samples="REF"), make_guide(samples="ALT")]
    grouped = group_guides_position(guides, debug=False)
    assert isinstance(grouped, dict)
    for k, v in grouped.items():
        assert isinstance(v, tuple)


def test_cfdon_score(monkeypatch):
    prepare_package()
    guides = [make_guide(samples="REF"), make_guide(samples="ALT")]
    monkeypatch.setattr("crisprhawk.scoring.cfdon", lambda ref, gs, debug: [0.9, 0.8])
    scored = cfdon_score(guides, verbosity=0, debug=False)
    assert all(g.cfdon_score in ("0.9", "0.8") for g in scored)


def test_elevationon_score(monkeypatch):
    prepare_package()
    guides = [make_guide(samples="REF"), make_guide(samples="ALT")]
    monkeypatch.setattr("crisprhawk.scoring.elevationon", lambda groups: guides)
    scored = elevationon_score(guides, verbosity=0, debug=False)
    assert scored == guides


def test_outofframe_score(monkeypatch):
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.scoring.ooframe_score", lambda guides, idx: [2])
    scored = outofframe_score(guides, guidelen=23, right=True, verbosity=0, debug=False)
    assert scored[0].ooframe_score == "2"