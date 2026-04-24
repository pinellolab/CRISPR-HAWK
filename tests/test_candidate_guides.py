import pytest
import pandas as pd
import numpy as np
import os
from crisprhawk import candidate_guides as cgmod
from crisprhawk.coordinate import Coordinate
from crisprhawk.pam import SPCAS9, XCAS9


class DummyRegion:
    def __init__(self, contig="chr1", start=100, end=200):
        self.coordinates = Coordinate(contig, start, end, 0)

    def contains(self, coord):
        return (
            self.coordinates.contig == coord.contig
            and self.coordinates.start <= coord.start <= self.coordinates.stop
        )


class DummyPAM:
    def __init__(self, name="NGG", cas_system=SPCAS9):
        self.name = name
        self.cas_system = cas_system

    def __str__(self):
        return self.name


@pytest.fixture
def dummy_cg():
    return cgmod.CandidateGuide("chr1:120:+", 20, False)


@pytest.fixture
def dummy_report_df():
    # Construct DF using the indices from REPORTCOLS
    cols = cgmod.REPORTCOLS
    data = {
        cols[1]: [120, 120, 150],  # Positions
        cols[10]: [0.95, 0.80, 0.75],  # CFD scores
        cols[13]: ["alt", "ref", "alt"],  # Guide type
        cols[14]: ["S1,S2", "S1", "S3"],  # Samples
    }
    return pd.DataFrame(data)


def test_candidate_guide_parsing():
    cg = cgmod.CandidateGuide("chr1:500:-", 20, False)
    assert cg.contig == "chr1"
    assert cg.position == 500
    assert cg.strand == "-"


def test_subset_region_report(dummy_report_df, dummy_cg):
    sub = cgmod._subset_region_report(dummy_report_df, dummy_cg, False)
    assert len(sub) == 2
    assert all(sub[cgmod.REPORTCOLS[1]] == 120)


def test_retrieve_scores_samples(dummy_report_df):
    # Only 'alt' rows should be retrieved
    scores, samples = cgmod._retrieve_scores_samples(dummy_report_df)
    assert len(scores) == 2
    assert samples == [2, 1]  # "S1,S2" -> 2, "S3" -> 1 (ref is skipped)


def test_subset_reports_integration(tmp_path, dummy_report_df, dummy_cg):
    # Setup files
    report_path = tmp_path / "region_report.tsv"
    dummy_report_df.to_csv(report_path, sep="\t", index=False)

    region = DummyRegion()
    region_reports = {region.coordinates: str(report_path)}
    pam = DummyPAM()

    result = cgmod.subset_reports([dummy_cg], region_reports, pam, 20, str(tmp_path), False)  # type: ignore

    assert dummy_cg in result
    assert os.path.exists(result[dummy_cg])
    # Check filename contains PAM name
    assert "NGG" in result[dummy_cg]


def test_prepare_data_dotplot(tmp_path, dummy_report_df, dummy_cg):
    report_path = tmp_path / "cg_report.tsv"
    dummy_report_df.to_csv(report_path, sep="\t", index=False)

    cg_reports = {dummy_cg: str(report_path)}
    data = cgmod._prepare_data_dotplot(cg_reports)

    label = str(dummy_cg)
    assert label in data
    scores, samples = data[label]
    assert len(scores) == 2  # Only alts
