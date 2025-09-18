from crisprhawk.annotate import (
    reverse_guides,
    azimuth_score,
    rs3_score,
    deepcpf1_score,
    group_guides_position,
    cfdon_score,
    elevationon_score,
    polish_variants_annotation,
    annotate_variants,
    annotate_variants_afs,
    gc_content,
    outofframe_score,
)
from crisprhawk.guide import Guide
from crisprhawk.region import Region
from crisprhawk.pam import PAM, SPCAS9, CPF1
from crisprhawk.config_crispritz import CrispritzConfig
from crisprhawk.utils import prepare_package

import numpy as np
import pytest

class DummyBedAnnotation:
    def __init__(self, *args, **kwargs):
        pass
    def fetch_features(self, contig, start, stop):
        return [f"{contig} {start} {stop} feature1 feature2 feature3 feature4 feature5 feature6 gene1;gene_name=GENE1;"]
    
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
        hapid="hap1"
    )
    defaults.update(kwargs)
    return Guide(**defaults)

def test_reverse_guides():
    guides = [make_guide(direction=1), make_guide(direction=0)]
    reversed_guides = reverse_guides(guides, verbosity=0)
    assert len(reversed_guides) == 2

def test_azimuth_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotate.azimuth", lambda seqs: [0.5])
    scored = azimuth_score(guides, verbosity=0, debug=False)
    assert scored[0].azimuth_score == "0.5"

def test_rs3_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotate.rs3", lambda seqs: [0.7])
    scored = rs3_score(guides, verbosity=0, debug=False)
    assert scored[0].rs3_score == "0.7"

def test_deepcpf1_score(monkeypatch):
    prepare_package()
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotate.deepcpf1", lambda seqs: [0.8])
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
    monkeypatch.setattr("crisprhawk.annotate.cfdon", lambda ref, gs, debug: [0.9, 0.8])
    scored = cfdon_score(guides, verbosity=0, debug=False)
    assert all(g.cfdon_score in ("0.9", "0.8") for g in scored)

def test_elevationon_score(monkeypatch):
    prepare_package()
    guides = [make_guide(samples="REF"), make_guide(samples="ALT")]
    monkeypatch.setattr("crisprhawk.annotate.elevationon", lambda groups: guides)
    scored = elevationon_score(guides, verbosity=0, debug=False)
    assert scored == guides

def test_polish_variants_annotation():
    guide = make_guide()
    variants = {"chr1-10-A/T", "chr1-12-G/C"}
    result = polish_variants_annotation(guide, variants)
    assert isinstance(result, set)

def test_annotate_variants_afs():
    guides = [make_guide(variants="chr1-10-A/T,chr1-12-G/C", afs={"chr1-10-A/T": 0.1, "chr1-12-G/C": 0.2})]
    annotated = annotate_variants_afs(guides, verbosity=0)
    assert isinstance(annotated, list)
    assert annotated[0].afs_str == "0.1,0.2"

def test_gc_content(monkeypatch):
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotate.gc_fraction", lambda seq: 0.5)
    gc_guides = gc_content(guides, verbosity=0, debug=False)
    assert gc_guides[0].gc == "0.5"

def test_outofframe_score(monkeypatch):
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotate.ooframe_score", lambda guides, idx: [2])
    scored = outofframe_score(guides, guidelen=23, right=True, verbosity=0, debug=False)
    assert scored[0].ooframe_score == "2"
