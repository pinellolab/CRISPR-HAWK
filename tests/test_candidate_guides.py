from crisprhawk.candidate_guides import (
    CandidateGuide,
    _parse_candidate_coord,
    initialize_candidate_guides,
    initialize_region_reports,
    _subset_region_report,
    _retrieve_scores_samples,
    _prepare_scatter_data,
    _extract_variant_labels,
)
from crisprhawk.coordinate import Coordinate

import pandas as pd
import numpy as np

import pytest
import os


class DummyRegion:
    def __init__(self, coordinates):
        self.coordinates = coordinates

    def contains(self, coord):
        return self.coordinates == coord


def test_candidate_guide_str_and_properties():
    cg = CandidateGuide("chr1:100", 20, False)
    assert str(cg) == "chr1:100"
    assert cg.contig == "chr1"
    assert cg.position == 100
    assert isinstance(cg.coordinate, Coordinate)


def test_parse_candidate_coord_valid():
    coord = _parse_candidate_coord("chr2:200", 30, False)
    assert coord.contig == "chr2"
    assert coord.start == 200
    assert coord.stop == 230


def test_initialize_candidate_guides():
    guides = ["chr1:100", "chr2:200"]
    result = initialize_candidate_guides(guides, 20, False)
    assert len(result) == 2
    assert all(isinstance(cg, CandidateGuide) for cg in result)


def test_initialize_region_reports():
    coord1 = Coordinate("chr1", 100, 120, 0)
    coord2 = Coordinate("chr2", 200, 220, 0)
    region1 = DummyRegion(coord1)
    region2 = DummyRegion(coord2)
    report_fnames = {region1: "file1.tsv", region2: "file2.tsv"}
    result = initialize_region_reports(report_fnames)
    assert result[coord1] == "file1.tsv"
    assert result[coord2] == "file2.tsv"


def test_subset_region_report(tmp_path):
    # Setup
    from crisprhawk.reports import REPORTCOLS

    cg = CandidateGuide("chr1:100", 20, False)
    data = {
        REPORTCOLS[1]: [100, 101, 102],
        REPORTCOLS[10]: [0.9, 0.8, 0.7],
        REPORTCOLS[13]: ["alt", "alt", "alt"],
        REPORTCOLS[14]: ["A,B", "A", "B"],
        REPORTCOLS[15]: ["var1", "var2", "var3"],
        REPORTCOLS[20]: [0.5, 0.6, 0.7],
    }
    df = pd.DataFrame(data)
    result = _subset_region_report(df, cg, False)
    assert not result.empty
    assert all(result[REPORTCOLS[1]] == 100)


def test_retrieve_scores_samples():
    from crisprhawk.reports import REPORTCOLS

    data = {
        REPORTCOLS[10]: ["0.5", "0.7", "0.9"],
        REPORTCOLS[13]: ["alt", "ref", "alt"],
        REPORTCOLS[14]: ["A,B", "A", "B,C,D"],
    }
    df = pd.DataFrame(data)
    scores, numsamples = _retrieve_scores_samples(df)
    assert scores == [0.5, 0.9]
    assert numsamples == [2, 3]


def test_prepare_scatter_data():
    from crisprhawk.reports import REPORTCOLS

    data = {
        REPORTCOLS[14]: ["A,B", "A", "B,C,D"],
    }
    df = pd.DataFrame(data)
    result = _prepare_scatter_data(df)
    assert "n_samples" in result.columns
    assert list(result["n_samples"]) == [2, 1, 3]


def test_extract_variant_labels():
    from crisprhawk.reports import REPORTCOLS

    data = {
        REPORTCOLS[15]: ["var1", float("nan"), "var3"],
    }
    df = pd.DataFrame(data)
    labels = _extract_variant_labels(df)
    assert labels == ["var1", "REF", "var3"]
