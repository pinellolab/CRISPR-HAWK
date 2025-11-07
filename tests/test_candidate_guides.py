import pytest
import pandas as pd
import numpy as np
import os
from types import SimpleNamespace

import crisprhawk.candidate_guides as cgmod
from crisprhawk.coordinate import Coordinate


class DummyRegion:
    def __init__(self, contig="chr1", start=100, end=120):
        self.coordinates = Coordinate(contig, start, end, 0)

    def contains(self, coord):
        return True


class DummyPAM:
    def __init__(self, name="NGG", cas_system="spcas9"):
        self.name = name
        self.cas_system = cas_system

    def __str__(self):
        return self.name


@pytest.fixture
def dummy_candidate_guide():
    return cgmod.CandidateGuide("chr1:100:+", 20, False)


@pytest.fixture
def dummy_report_df():
    # Minimal DataFrame with required columns for candidate guide analysis
    data = {
        cgmod.REPORTCOLS[1]: [100, 100, 101],
        cgmod.REPORTCOLS[10]: [0.9, 0.8, 0.7],
        cgmod.REPORTCOLS[13]: ["alt", "ref", "alt"],
        cgmod.REPORTCOLS[14]: ["sample1,sample2", "sample1", "sample3"],
        cgmod.REPORTCOLS[15]: ["var1", np.nan, "var2"],
        cgmod.REPORTCOLS[20]: [0.5, 0.6, 0.7],
    }
    return pd.DataFrame(data)


def test_initialize_candidate_guides(dummy_candidate_guide):
    guides = cgmod.initialize_candidate_guides(["chr1:100:+", "chr1:200:-"], 20, False)
    assert all(isinstance(g, cgmod.CandidateGuide) for g in guides)


def test_initialize_region_reports():
    region1 = DummyRegion()
    region2 = DummyRegion("chr2", 200, 220)
    reports = {region1: "file1.tsv", region2: "file2.tsv"}
    mapping = cgmod.initialize_region_reports(reports)
    assert all(isinstance(k, Coordinate) for k in mapping.keys())
    assert set(mapping.values()) == {"file1.tsv", "file2.tsv"}


def test_subset_region_report_found(dummy_report_df, dummy_candidate_guide):
    # Should return a DataFrame with only rows matching candidate guide position
    result = cgmod._subset_region_report(dummy_report_df, dummy_candidate_guide, False)
    assert not result.empty
    assert all(result[cgmod.REPORTCOLS[1]] == 100)


def test_subset_region_report_not_found(dummy_report_df):
    # Should call exception_handler if not found
    dummy_cg = cgmod.CandidateGuide("chr1:999:+", 20, False)
    called = {}

    def fake_exception_handler(*args, **kwargs):
        called["called"] = True
        raise Exception("Handled")

    cgmod.exception_handler = fake_exception_handler
    with pytest.raises(Exception):
        cgmod._subset_region_report(dummy_report_df, dummy_cg, False)
    assert called.get("called", False)


def test_store_region_report_subset(tmp_path, dummy_report_df, dummy_candidate_guide):
    fname = cgmod._store_region_report_subset(
        dummy_report_df, dummy_candidate_guide, DummyPAM(), 20, str(tmp_path)
    )
    assert os.path.exists(fname)
    df = pd.read_csv(fname, sep="\t")
    assert not df.empty


def test_subset_reports(tmp_path, dummy_report_df, dummy_candidate_guide):
    # Save dummy report to file
    report_path = os.path.join(tmp_path, "region.tsv")
    dummy_report_df.to_csv(report_path, sep="\t", index=False)
    region = DummyRegion()
    region_reports = {region.coordinates: report_path}
    result = cgmod.subset_reports(
        [dummy_candidate_guide], region_reports, DummyPAM(), 20, str(tmp_path), False
    )
    assert dummy_candidate_guide in result
    assert os.path.exists(result[dummy_candidate_guide])


def test_retrieve_scores_samples(dummy_report_df):
    scores, numsamples = cgmod._retrieve_scores_samples(dummy_report_df)
    assert isinstance(scores, list)
    assert isinstance(numsamples, list)
    assert all(isinstance(s, float) for s in scores)
    assert all(isinstance(n, int) for n in numsamples)


def test_prepare_data_dotplot(tmp_path, dummy_report_df, dummy_candidate_guide):
    # Save dummy report to file
    report_path = os.path.join(tmp_path, "cg.tsv")
    dummy_report_df.to_csv(report_path, sep="\t", index=False)
    cg_reports = {dummy_candidate_guide: report_path}
    data = cgmod._prepare_data_dotplot(cg_reports)
    assert isinstance(data, dict)
    assert str(dummy_candidate_guide) in data


def test_setup_figure():
    fig, ax = cgmod._setup_figure()
    assert fig is not None
    assert ax is not None


def test_color_labels_dotplot():
    labels = ["cg1", "cg2", "cg3"]
    colorlab = cgmod._color_labels_dotplot(labels)
    assert set(colorlab.keys()) == set(labels)


def test_gradient_colormap():
    cmap = cgmod._gradient_colormap()
    assert hasattr(cmap, "__call__")


def test_background_grid_dotplot():
    arr = cgmod._background_grid_dotplot(["cg1", "cg2"])
    assert isinstance(arr, np.ndarray)
    assert arr.shape[0] > 0


def test_draw_background_gradient():
    labels = ["cg1", "cg2"]
    cmap = cgmod._gradient_colormap()
    im = cgmod._draw_background_gradient(labels, cmap)
    assert hasattr(im, "get_array")


def test_plot_candidate_guide_data_runs():
    # Just check that it runs
    cgmod._plot_candidate_guide_data("cg1", [0.8, 0.9], [1, 2], 0, (0.5, 0.5, 0.5))


def test_plot_all_data_points_runs():
    dotplot_data = {"cg1": ([0.8, 0.9], [1, 2]), "cg2": ([0.7], [3])}
    labels = ["cg1", "cg2"]
    colorlabs = {k: (0.5, 0.5, 0.5) for k in labels}
    cgmod._plot_all_data_points(dotplot_data, labels, colorlabs)


def test_add_reference_line_runs():
    cgmod._add_reference_line()


def test_draw_legend_dotplot():
    handles = cgmod._draw_legend_dotplot()
    assert isinstance(handles, list)
    assert any(hasattr(h, "get_label") for h in handles)


def test_configure_legend_runs():
    cgmod._configure_legend()


def test_configure_axes_style_runs():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    cgmod._configure_axes_style(["cg1", "cg2"], ax)


def test_add_colorbar_runs():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    cmap = cgmod._gradient_colormap()
    im = cgmod._draw_background_gradient(["cg1", "cg2"], cmap)
    cgmod._add_colorbar(im, ax)


def test_save_figure_runs(tmp_path):
    cgmod._save_figure(str(tmp_path))


def test_draw_dotplot_runs(tmp_path):
    dotplot_data = {"cg1": ([0.8, 0.9], [1, 2]), "cg2": ([0.7], [3])}
    cgmod._draw_dotplot(dotplot_data, str(tmp_path))


def test_candidate_guides_dotplot_runs(tmp_path):
    dotplot_data = {"cg1": ([0.8, 0.9], [1, 2]), "cg2": ([0.7], [3])}
    cgmod.candidate_guides_dotplot(dotplot_data, str(tmp_path))


def test_prepare_scatter_data(dummy_report_df):
    df = cgmod._prepare_scatter_data(dummy_report_df)
    assert "n_samples" in df.columns


def test_extract_variant_labels(dummy_report_df):
    labels = cgmod._extract_variant_labels(dummy_report_df)
    assert isinstance(labels, list)
    assert len(labels) == len(dummy_report_df)


def test_create_diagonal_gradient():
    arr, cmap = cgmod._create_diagonal_gradient()
    assert isinstance(arr, np.ndarray)
    assert hasattr(cmap, "__call__")


def test_draw_background_scatter_runs():
    arr, cmap = cgmod._create_diagonal_gradient()
    im = cgmod._draw_background_scatter(arr, cmap)
    assert hasattr(im, "get_array")


def test_plot_scatter_points_runs(dummy_report_df):
    colorlab = {
        "var1": (0.5, 0.5, 0.5),
        "var2": (0.2, 0.2, 0.2),
        "REF": (0.1, 0.1, 0.1),
    }
    dummy_report_df = cgmod._prepare_scatter_data(dummy_report_df)
    cgmod._plot_scatter_points(dummy_report_df, colorlab)


def test_create_sample_size_legend_runs():
    handles = cgmod._create_sample_size_legend([1, 2, 3])
    assert isinstance(handles, list)


def test_configure_scatter_axes_runs(dummy_candidate_guide):
    cgmod._configure_scatter_axes(dummy_candidate_guide)


def test_add_scatter_colorbar_runs():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    arr, cmap = cgmod._create_diagonal_gradient()
    im = cgmod._draw_background_scatter(arr, cmap)
    cgmod._add_scatter_colorbar(im, fig, ax)


def test_add_scatter_legends_runs(dummy_report_df):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    colorlab = {
        "var1": (0.5, 0.5, 0.5),
        "var2": (0.2, 0.2, 0.2),
        "REF": (0.1, 0.1, 0.1),
    }
    dummy_report_df = cgmod._prepare_scatter_data(dummy_report_df)
    variant_handles = cgmod._plot_scatter_points(dummy_report_df, colorlab)
    size_handles = cgmod._create_sample_size_legend([1, 2, 3])
    cgmod._add_scatter_legends(variant_handles, size_handles)


def test_save_scatter_figure_runs(tmp_path, dummy_candidate_guide):
    cgmod._save_scatter_figure(dummy_candidate_guide, str(tmp_path))


def test_draw_scatter_runs(tmp_path, dummy_report_df, dummy_candidate_guide):
    cgmod._draw_scatter(dummy_report_df, dummy_candidate_guide, str(tmp_path))


def test_candidate_guides_scatter_runs(
    tmp_path, dummy_report_df, dummy_candidate_guide
):
    # Save dummy report to file
    report_path = os.path.join(tmp_path, "cg.tsv")
    dummy_report_df.to_csv(report_path, sep="\t", index=False)
    cg_reports = {dummy_candidate_guide: report_path}
    cgmod.candidate_guides_scatter(cg_reports, str(tmp_path))


def test_draw_candidate_guides_plots_runs(
    tmp_path, dummy_report_df, dummy_candidate_guide
):
    # Save dummy report to file
    report_path = os.path.join(tmp_path, "cg.tsv")
    dummy_report_df.to_csv(report_path, sep="\t", index=False)
    cg_reports = {dummy_candidate_guide: report_path}
    cgmod.draw_candidate_guides_plots(cg_reports, True, str(tmp_path))


def test_candidate_guides_analysis_runs(tmp_path, dummy_report_df):
    # Save dummy report to file
    report_path = os.path.join(tmp_path, "region.tsv")
    dummy_report_df.to_csv(report_path, sep="\t", index=False)
    region = DummyRegion()
    reports = {region: report_path}
    pam = DummyPAM()
    args = SimpleNamespace(
        candidate_guides=["chr1:100:+"],
        guidelen=20,
        debug=False,
        outdir=str(tmp_path),
        estimate_offtargets=True,
        verbosity=1,
    )
    cgmod.candidate_guides_analysis(reports, pam, args)
