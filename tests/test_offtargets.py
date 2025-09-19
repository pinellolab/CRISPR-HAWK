from crisprhawk.offtargets import (
    _write_guides_file,
    _write_pam_file,
    _prepare_input_data,
    _format_targets_prefix,
    _calculate_offtargets_map,
    _calculate_global_cfd,
    annotate_guides_offtargets,
)
from crisprhawk.guide import Guide
from crisprhawk.pam import PAM
from crisprhawk.region import Region

import tempfile
import pytest
import os

class DummyConfig:
    outdir = ".crispritz_targets"
    conda = "conda"
    env_name = "crispritz"

def test_write_guides_file_creates_file():
    guides = [Guide(0, 120, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGACTGACTGACTGTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 20, 3, 0, "s1", "NA", {}, False, False, "hap1")]
    pam = PAM("NGG", False, False)
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = _write_guides_file(guides, pam, tmpdir, verbosity=0, debug=False)
        assert os.path.exists(fname)
        with open(fname) as f:
            lines = f.readlines()
        assert len(lines) == 1
        assert "ACTGACTGACTGACTGACTG" in lines[0].strip()

def test_write_pam_file_creates_file():
    pam = PAM("NGG", False, False)
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = _write_pam_file(pam, guidelen=20, right=False, crispritz_dir=tmpdir, verbosity=0, debug=False)
        assert os.path.exists(fname)
        with open(fname) as f:
            lines = f.readlines()
        assert len(lines) == 1
        assert "NGG" in lines[0]

def test_prepare_input_data_creates_files():
    guides = [Guide(0, 120, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGACTGACTGACTGTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 20, 3, 0, "s1", "NA", {}, False, False, "hap1")]
    pam = PAM("NGG", False, False)
    config = DummyConfig()
    with tempfile.TemporaryDirectory() as tmpdir:
        guides_fname, pam_fname = _prepare_input_data(config, guides, pam, tmpdir, verbosity=0, debug=False)
        assert os.path.exists(guides_fname)
        assert os.path.exists(pam_fname)

def test_calculate_offtargets_map_returns_dict():
    class DummyOfftarget:
        def __init__(self, grna_):
            self.grna_ = grna_
    guides = [Guide(0, 120, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGACTGACTGACTGTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 20, 3, 0, "s1", "NA", {}, False, False, "hap1")]
    offtargets = [DummyOfftarget("ACTGACTGACTGACTGACTG")]
    otmap = _calculate_offtargets_map(offtargets, guides)
    assert isinstance(otmap, dict)
    assert guides[0].guide.upper() in otmap

def test_calculate_global_cfd_returns_float():
    class DummyOfftarget:
        def __init__(self, cfd):
            self.cfd = cfd
    offtargets = [DummyOfftarget("0.5"), DummyOfftarget("NA")]
    score = _calculate_global_cfd(offtargets, verbosity=0)
    assert isinstance(score, float)

def test_annotate_guides_offtargets_sets_attributes():
    class DummyOfftarget:
        def __init__(self, grna_, cfd):
            self.grna_ = grna_
            self.cfd = cfd
    guides = guides = [Guide(0, 120, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACTGACTGACTGACTGACTGTGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", 20, 3, 0, "s1", "NA", {}, False, False, "hap1")]
    offtargets = [DummyOfftarget("ACTGACTGACTGACTGACTG", "0.5")]
    guides = annotate_guides_offtargets(offtargets, guides, verbosity=0)
    assert hasattr(guides[0], "offtargets")
    assert hasattr(guides[0], "cfd")