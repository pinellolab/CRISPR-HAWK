""" """

from .crisprhawk_error import CrisprHawkCandidateGuideError
from .exception_handlers import exception_handler
from .coordinate import Coordinate
from .utils import warning, print_verbosity, CANDIDATEGUIDESREPORTPREFIX, VERBOSITYLVL
from .graphical_reports import FIGURE_SIZE
from .reports import REPORTCOLS
from .region import Region
from .pam import PAM, SPCAS9, XCAS9

from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.pyplot import Axes
from matplotlib.figure import Figure
from matplotlib.cm import ScalarMappable
from typing import List, Dict, Tuple
from time import time

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np

import os

GRADIENTCOLORS = ["#B30000", "#FF4444", "#FFAAAA", "#FFFFFF"]
DOTPLOTLEGENDLABS = [
    "1 sample",
    "2-20 sample",
    "21-50 sample",
    "51-100 sample",
    "101-200 sample",
    "> 200 sample",
]


class CandidateGuide:

    def __init__(self, coordinate: str, guidelen: int, debug: bool) -> None:
        self._debug = debug  # store debug flag
        self._coordinate = _parse_candidate_coord(coordinate, guidelen, self._debug)

    def __str__(self) -> str:
        return f"{self.contig}:{self.position}"

    @property
    def coordinate(self) -> Coordinate:
        return self._coordinate

    @property
    def contig(self) -> str:
        return self._coordinate.contig

    @property
    def position(self) -> int:
        return self._coordinate.start


def _parse_candidate_coord(coordinate: str, guidelen: int, debug: bool) -> Coordinate:
    contig, position = coordinate.split(":")  # retrieve contig and position
    try:
        return Coordinate(contig, int(position), int(position) + guidelen, 0)
    except Exception as e:
        exception_handler(
            CrisprHawkCandidateGuideError,
            f"Forbidden candidate guide coordinate ({coordinate})",
            os.EX_DATAERR,
            debug,
            e,
        )


def initialize_candidate_guides(
    candidate_guides: List[str], guidelen: int, debug: bool
) -> List[CandidateGuide]:
    """Creates a list of CandidateGuide objects from provided guide coordinates.

    This function initializes CandidateGuide instances for each guide coordinate
    string in the input list.

    Args:
        candidate_guides: List of guide coordinate strings.
        guidelen: Length of each guide.
        debug: Flag to enable debug mode.

    Returns:
        List of initialized CandidateGuide objects.
    """
    return [CandidateGuide(g, guidelen, debug) for g in candidate_guides]


def initialize_region_reports(
    report_fnames: Dict[Region, str],
) -> Dict[Coordinate, str]:
    """Creates a mapping from region coordinates to report filenames.

    This function generates a dictionary that links each region's coordinates to
    its corresponding report file.

    Args:
        report_fnames: Dictionary mapping Region objects to report filenames.

    Returns:
        Dictionary mapping Coordinate objects to report filenames.
    """
    return {region.coordinates: fname for region, fname in report_fnames.items()}


def _subset_region_report(
    report: pd.DataFrame, cg: CandidateGuide, debug: bool
) -> pd.DataFrame:
    """Extracts the subset of a region report corresponding to a candidate guide.

    This function filters the region report to include only entries matching the
    candidate guide's position.

    Args:
        report: DataFrame containing the region report.
        cg: CandidateGuide object for which to extract the subset.
        debug: Flag to enable debug mode.

    Returns:
        DataFrame containing only the rows relevant to the candidate guide.

    Raises:
        CrisprHawkCandidateGuideError: If the candidate guide is not found in
            the report.
    """
    report_sub = report[
        report[REPORTCOLS[1]] == cg.position
    ]  # only alternative to candidate
    if report_sub.empty:  # empty?
        exception_handler(
            CrisprHawkCandidateGuideError,
            f"Candidate guide {cg} not found. Is the candidate guide correct?",
            os.EX_DATAERR,
            debug,
        )
    report_sub.reset_index(drop=True, inplace=True)  # reset index
    return report_sub


def _store_region_report_subset(
    report_sub: pd.DataFrame, cg: CandidateGuide, pam: PAM, guidelen: int, outdir: str
) -> str:
    """Stores a subset of a region report for a candidate guide to a file.

    This function saves the provided DataFrame as a TSV file named according to
    the candidate guide and parameters.

    Args:
        report_sub: DataFrame containing the subset of the region report.
        cg: CandidateGuide object for which the report is generated.
        pam: PAM object specifying the protospacer adjacent motif.
        guidelen: Length of the guide.
        outdir: Output directory for the report file.

    Returns:
        The filename of the stored subset report.
    """
    report_sub_fname = os.path.join(
        outdir,
        f"{CANDIDATEGUIDESREPORTPREFIX}__{cg.contig}_{cg.position}_{pam}_{guidelen}.tsv",
    )
    report_sub.to_csv(
        report_sub_fname, sep="\t", index=False, na_rep="NA"
    )  # store sub report
    return report_sub_fname


def subset_reports(
    candidate_guides: List[CandidateGuide],
    region_reports: Dict[Coordinate, str],
    pam: PAM,
    guidelen: int,
    outdir: str,
    debug: bool,
) -> Dict[CandidateGuide, str]:
    """Generates and stores subset reports for candidate guides within specified regions.

    This function processes region reports to extract and save data relevant to
    each candidate guide. It returns a mapping from each CandidateGuide to its
    corresponding subset report file.

    Args:
        candidate_guides: List of CandidateGuide objects to analyze.
        region_reports: Dictionary mapping region coordinates to report filenames.
        pam: PAM object specifying the protospacer adjacent motif.
        guidelen: Length of each guide.
        outdir: Output directory for storing subset reports.
        debug: Flag to enable debug mode.

    Returns:
        Dictionary mapping CandidateGuide objects to their subset report filenames.
    """
    cg_reports = {}  # candidate guides sub reports
    for region, report_fname in region_reports.items():
        report = pd.read_csv(report_fname, sep="\t")  # load region grnas report
        for cg in candidate_guides:  # retrieve candidate guides in current region
            if region.contains(cg.coordinate):
                report_sub = _subset_region_report(report, cg, debug)
                if not report_sub.empty:
                    cg_reports[cg] = _store_region_report_subset(
                        report_sub, cg, pam, guidelen, outdir
                    )
    return cg_reports


def _retrieve_scores_samples(report: pd.DataFrame) -> Tuple[List[float], List[int]]:
    """Extracts scores and sample counts from a candidate guide report.

    This function filters out reference guides and retrieves CFD scores and the
    number of samples for each guide.

    Args:
        report: DataFrame containing the candidate guide report.

    Returns:
        Tuple containing a list of CFD scores and a list of sample counts.
    """
    rdf = report[report[REPORTCOLS[14]] != "ref"]  # keep only alternative grnas
    scores = list(map(float, rdf[REPORTCOLS[10]].tolist()))  # cfdon scores
    # count number of samples per grna
    numsamples = rdf[REPORTCOLS[15]].str.split(",").apply(len)
    numsamples = numsamples.fillna(1)  # fill potential na values
    return scores, list(numsamples)


def _prepare_data_dotplot(
    cg_reports: Dict[CandidateGuide, str], debug: bool
) -> Dict[str, Tuple[List[float], List[int]]]:
    """Prepares data for dot plot visualization of candidate guide variant effects.

    This function reads each candidate guide's report and extracts scores and sample
    counts for plotting.

    Args:
        cg_reports: Dictionary mapping CandidateGuide objects to their subset
            report filenames.
        debug: Flag to enable debug mode.

    Returns:
        Dictionary mapping candidate guide string representations to tuples of
            scores and sample counts.
    """
    dotplot_data = {}  # dotplot data dictionary
    for cg, report in cg_reports.items():
        rdf = pd.read_csv(report, sep="\t")  # read candidate guide report
        # retrieve scores and number of samples
        dotplot_data[f"{cg}"] = _retrieve_scores_samples(rdf)
    return dotplot_data


def _setup_figure() -> Tuple[Figure, Axes]:
    """Sets up a matplotlib figure and axes for dot plot visualization.

    This function creates and returns a new figure and axes with a predefined
    size for plotting.

    Returns:
        Tuple containing the Figure and Axes objects.
    """
    return plt.subplots(1, 1, figsize=(25, 15))


def _color_labels_dotplot(labels: List[str]) -> Dict[str, Tuple]:
    """Assigns unique colors to each candidate guide label for dot plot visualization.

    This function generates a color mapping for each label using a categorical colormap.

    Args:
        labels: List of candidate guide string labels.

    Returns:
        Dictionary mapping each label to a color tuple.
    """
    cmap = cm.get_cmap("tab20", len(labels))  # compute colormap
    return {label: cmap(i) for i, label in enumerate(labels)}


def _gradient_colormap() -> LinearSegmentedColormap:
    """Creates a custom linear segmented colormap for dot plot backgrounds.

    This function generates a gradient colormap using predefined color stops.

    Returns:
        A LinearSegmentedColormap object for use in background gradients.
    """
    return LinearSegmentedColormap.from_list("custom_gradient", GRADIENTCOLORS, N=256)


def _background_grid_dotplot(labels: List[str]) -> np.ndarray:
    """Generates a vertical gradient grid for the dot plot background.

    This function creates a normalized 2D array representing the vertical gradient
    for the plot background.

    Args:
        labels: List of candidate guide string labels.

    Returns:
        A 2D numpy array with normalized values for the vertical gradient.
    """
    x = np.linspace(-0.5, len(labels) - 0.5, 256)
    y = np.linspace(-0.1, 1.1, 256)
    X, Y = np.meshgrid(x, y)
    # Normalize Y to create the vertical gradient
    return (Y - Y.min()) / (Y.max() - Y.min())


def _draw_background_gradient(
    labels: List[str], colorsgrads: LinearSegmentedColormap
) -> ScalarMappable:
    """Draws a background gradient for the dot plot visualization.

    This function creates a vertical gradient background for the dot plot using
    the provided colormap.

    Args:
        labels: List of candidate guide string labels.
        colorsgrads: LinearSegmentedColormap object for the gradient.

    Returns:
        ScalarMappable object representing the background image.
    """
    vdistance = _background_grid_dotplot(labels)
    return plt.imshow(
        vdistance,
        extent=[-0.5, len(labels) - 0.5, -0.1, 1.1],
        origin="lower",
        cmap=colorsgrads,
        aspect="auto",
        alpha=0.3,
        zorder=0,
    )


def _plot_candidate_guide_data(
    cg: str, scores: List[float], numsamples: List[int], x_position: int, color: Tuple
) -> None:
    """Plots data points for a single candidate guide on the dot plot.

    This function visualizes each score and sample count for the candidate guide
    at the specified x-axis position.

    Args:
        cg: String representation of the candidate guide.
        scores: List of CFD scores for the candidate guide.
        numsamples: List of sample counts corresponding to each score.
        x_position: X-axis position for plotting the candidate guide.
        color: Color tuple for the data points.

    Returns:
        None
    """
    # scale dot size by square root of sample count for better visual perception
    dotsizes = [150 * np.sqrt(ns) for ns in numsamples]
    for i, score in enumerate(scores):
        marker = "D" if numsamples[i] == 1 else "o"  # diamond for singletons
        plt.scatter(x_position, score, alpha=0.6, c=[color], s=dotsizes[i], marker=marker, edgecolors="white", linewidth=0.5, zorder=2)  # type: ignore


def _plot_all_data_points(
    dotplot_data: Dict[str, Tuple[List[float], List[int]]],
    labels: List[str],
    colorlabs: Dict[str, Tuple],
) -> None:
    """Plots all data points for candidate guides on the dot plot.

    This function iterates through the dotplot data and plots each candidate
    guide's scores and sample counts.

    Args:
        dotplot_data: Dictionary mapping candidate guide string representations
            to tuples of scores and sample counts.
        labels: List of candidate guide string labels.
        colorlabs: Dictionary mapping each label to a color tuple.

    Returns:
        None
    """
    labels_array = np.array(labels)
    for cg, data in dotplot_data.items():
        scores, numsamples = data
        x_position = int(np.where(labels_array == cg)[0][0])
        _plot_candidate_guide_data(cg, scores, numsamples, x_position, colorlabs[cg])


def _add_reference_line() -> None:
    """Adds a horizontal reference line at y=1 to the current plot.

    This function draws a gray horizontal line at y=1 to serve as a visual
    reference in the dot plot.

    Returns:
        None
    """
    plt.axhline(y=1, color="gray", linestyle="-", linewidth=3, alpha=0.6, zorder=1)


def _draw_legend_dotplot() -> List[Line2D]:
    """Creates legend handles for the dot plot visualization.

    This function generates custom legend handles representing sample size
    categories and the reference line.

    Returns:
        List of Line2D objects to be used as legend handles in the plot.
    """
    sample_counts = [1, 10, 35, 75, 150, 300]
    markers = ["D"] + ["o"] * 5
    scaled_sizes = [150 * np.sqrt(n) for n in sample_counts]  # scale legend dots size
    refhandle = Line2D(
        [0], [0], color="gray", linestyle="-", linewidth=3, label="Reference", alpha=0.6
    )
    handles_size = [
        Line2D(
            [0],
            [0],
            marker=marker,
            color="w",
            label=label,
            markerfacecolor="gray",
            markersize=np.sqrt(size),
            alpha=0.6,
        )
        for label, size, marker in zip(DOTPLOTLEGENDLABS, scaled_sizes, markers)
    ]
    return [refhandle] + handles_size


def _configure_legend() -> None:
    """Configures and displays the legend for the dot plot visualization.

    This function creates and adds a custom legend to the plot, including sample
    size markers and a reference line.

    Returns:
        None
    """
    handles = _draw_legend_dotplot()
    plt.legend(
        handles=handles,
        frameon=False,
        labelspacing=1.5,
        ncol=len(handles),
        loc="lower center",
        bbox_to_anchor=(0.5, -0.2),
        prop={"size": 14},
    )
    plt.subplots_adjust(bottom=0.15)


def _configure_axes_style(labels: List[str], ax: Axes) -> None:
    """Configures the style and appearance of the axes for the dot plot visualization.

    This function sets titles, labels, spines, ticks, limits, and grid style
    for the plot axes.

    Args:
        labels: List of candidate guide string labels.
        ax: Matplotlib Axes object to configure.

    Returns:
        None
    """
    # set titles and labels
    plt.title("Variant Effect on On-Targets in Candidate Guides", fontsize=18)
    plt.xlabel("Candidate guides", fontsize=14)
    plt.ylabel("Variant Effect on On-Targets (CFD)", fontsize=14)
    # remove top and right spines for cleaner look
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # configure ticks and limits
    plt.xticks(range(len(labels)), labels, fontsize=15, rotation=45, ha="right")
    plt.yticks(fontsize=15)
    plt.ylim(-0.1, 1.1)
    plt.xlim(-0.5, len(labels) - 0.5)
    plt.grid(True, alpha=0.3, linestyle="--", zorder=1, color="black")  # add grid


def _add_colorbar(im: ScalarMappable, ax: Axes) -> None:
    """Adds a colorbar to the plot to indicate guide penalty levels.

    This function creates and customizes a colorbar for the dot plot, labeling
    it to show the range of guide penalties.

    Args:
        im: ScalarMappable object representing the background image.
        ax: Matplotlib Axes object to which the colorbar will be added.

    Returns:
        None
    """
    cbar = plt.colorbar(im, ax=ax, fraction=0.015, pad=0.02, shrink=0.6)
    cbar.set_label(
        "Guide Penalty\n(Low → High)", rotation=270, labelpad=25, fontsize=12
    )
    cbar.set_ticks([])  # remove tick marks for cleaner appearance


def _save_figure(outdir: str) -> None:
    """Saves the current matplotlib figure to the specified output directory.

    This function writes the figure as a PNG file and closes the plot to free
    resources.

    Args:
        outdir: Output directory where the figure will be saved.

    Returns:
        None
    """
    filename = "candidate_guides_variant_effect_ontarget.png"
    output_path = os.path.join(outdir, filename)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, format="png", bbox_inches="tight")
    plt.close()


def _draw_dotplot(
    dotplot_data: Dict[str, Tuple[List[float], List[int]]], outdir: str
) -> None:
    """Draws and saves a dot plot for candidate guide variant effects.

    This function visualizes the variant effect data for candidate guides and
    saves the resulting plot to a file.

    Args:
        dotplot_data: Dictionary mapping candidate guide string representations
            to tuples of scores and sample counts.
        outdir: Output directory for saving the figure.

    Returns:
        None
    """
    # prepare data and color schemes
    labels = list(dotplot_data.keys())
    colorlabs = _color_labels_dotplot(labels)
    colorsgrads = _gradient_colormap()
    f, ax = _setup_figure()  # set up figure
    im = _draw_background_gradient(labels, colorsgrads)  # draw background gradient
    _plot_all_data_points(dotplot_data, labels, colorlabs)  # plot all data points
    _add_reference_line()  # add reference line
    # configure visual elements
    _configure_legend()
    _configure_axes_style(labels, ax)
    _add_colorbar(im, ax)
    _save_figure(outdir)  # save the figure


def candidate_guides_dotplot(
    cg_reports: Dict[CandidateGuide, str], outdir: str, debug: bool
) -> None:
    """Generates and saves a dot plot visualizing variant effects on candidate
    guides.

    This function creates a graphical report showing the distribution of variant
    effects for each candidate guide.

    Args:
        cg_reports: Dictionary mapping CandidateGuide objects to their subset
            report filenames.
        outdir: Output directory for saving the figure.
        debug: Flag to enable debug mode.

    Returns:
        None
    """
    outdir_gr = os.path.join(outdir, "figures")  # graphical report folder
    if not os.path.isdir(outdir_gr):  # create figures folder
        os.makedirs(outdir_gr)
    dotplot_data = _prepare_data_dotplot(cg_reports, debug)  # prepare data fro dotplot
    _draw_dotplot(dotplot_data, outdir_gr)  # draw variant effect dotplot


def candidate_guides_analysis(
    candidate_guides_coords: List[str],
    reports: Dict[Region, str],
    pam: PAM,
    guidelen: int,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    print_verbosity("Analyzing candidate guides", verbosity, VERBOSITYLVL[1])
    start = time()
    candidate_guides = initialize_candidate_guides(
        candidate_guides_coords, guidelen, debug
    )  # initialize canididate guide objects
    region_reports = initialize_region_reports(reports)  # create report map
    # subset reports for candidate guides
    cg_reports = subset_reports(
        candidate_guides, region_reports, pam, guidelen, outdir, debug
    )
    if pam.cas_system in [SPCAS9, XCAS9]:
        candidate_guides_dotplot(cg_reports, outdir, debug)
    print_verbosity(
        f"Candidate guides analysis completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
