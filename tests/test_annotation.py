from crisprhawk.annotation import (
    reverse_guides,
    polish_variants_annotation,
    annotate_variants_afs,
    gc_content,
)
from crisprhawk.guide import Guide

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

def test_reverse_guides():
    guides = [make_guide(direction=1), make_guide(direction=0)]
    reversed_guides = reverse_guides(guides, verbosity=0)
    assert len(reversed_guides) == 2

def test_polish_variants_annotation():
    guide = make_guide()
    variants = {"chr1-10-A/T", "chr1-12-G/C"}
    result = polish_variants_annotation(guide, variants)
    assert isinstance(result, set)


def test_annotate_variants_afs():
    guides = [
        make_guide(
            variants="chr1-10-A/T,chr1-12-G/C",
            afs={"chr1-10-A/T": 0.1, "chr1-12-G/C": 0.2},
        )
    ]
    annotated = annotate_variants_afs(guides, verbosity=0)
    assert isinstance(annotated, list)
    assert annotated[0].afs_str == "0.1,0.2"

def test_gc_content(monkeypatch):
    guides = [make_guide()]
    monkeypatch.setattr("crisprhawk.annotation.gc_fraction", lambda seq: 0.5)
    gc_guides = gc_content(guides, verbosity=0, debug=False)
    assert gc_guides[0].gc == "0.5"
