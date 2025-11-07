"""Provides functions for generating graphical reports for CRISPR guide analysis.

This module includes utilities for creating pie charts and delta dot plots to
visualize guide type distributions and score variations across genomic regions.
It supports processing guide data, ranking guides, and saving publication-ready
figures for downstream analysis.
"""

from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .crisprhawk_error import CrisprHawkGraphicalReportsError
from .region_constructor import PADDING
from .exception_handlers import exception_handler
from .utils import VERBOSITYLVL, is_lowercase, create_folder, print_verbosity, warning
from .candidate_guides import CandidateGuide, initialize_candidate_guides
from .reports import REPORTCOLS
from .region import Region


from typing import Dict, List, Tuple, Any, Set, Optional
from matplotlib.lines import Line2D
from matplotlib.pyplot import Axes
from matplotlib.figure import Figure
from time import time

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import random
import os


# constants

# guide types (pie chart)
GUIDETYPES = {
    0: "Reference Guides",
    1: "Spacer+PAM Alternative Guides",
    2: "Spacer Alternative Guides",
    3: "PAM Alternative Guides",
}

# supported scores (dotplot)
SCORES = [
    "score_azimuth",
    "score_rs3",
    "score_deepcpf1",
    "score_cfdon",
    "score_elevationon",
]

# figures size
FIGURE_SIZE = (25, 20)

# figures dpi
DPI = 300

# base colormap
BASE_CMAPS = [
    "Purples",
    "Blues",
    "Greens",
    "Oranges",
    "Reds",
    "PuRd",
    "RdPu",
    "BuPu",
    "GnBu",
    "PuBuGn",
    "BuGn",
    "Spectral",
    "coolwarm",
]

# sample dots size multiplier
SAMPLE_SIZE_MULTIPLIER = 150

# representative samples
REPRESENTATIVE_SAMPLES = 30

# delta dotplot colnames
DELTACOLS = ["delta", "abs_delta"]

# legend configuration
LEGEND_LABELS = ["1", "2-20", "21-50", "51-100", "101-200", ">200"]
LEGEND_SAMPLE_COUNTS = [1, 10, 35, 75, 150, 300]
LEGEND_MARKERS = ["D", "o", "o", "o", "o", "o"]


def create_figures_dir(outdir: str) -> str:
    """Creates a 'figures' subdirectory in the specified output directory if it
    does not exist.

    Returns the path to the figures directory for saving graphical report images.

    Args:
        outdir: The base output directory where the figures folder will be created.

    Returns:
        str: The path to the created or existing figures directory.
    """
    # create figures folder in output directory
    outdir_gr = os.path.join(outdir, "figures")
    if not os.path.isdir(outdir_gr):
        outdir_gr = create_folder(os.path.join(outdir, "figures"))
    return outdir_gr


def format_region_prefix(region: Region) -> str:
    """Formats a Region object into a string for use in plot naming and reporting.

    The returned string includes the contig, padded start, and padded stop positions.

    Args:
        region: A Region object containing contig, start, and stop attributes.

    Returns:
        str: A formatted string representing the region.
    """
    return f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}"


def _compute_extended_guide_id(
    chrom: str, start: int, stop: int, strand: str, sgrna: str, pam: str
) -> str:
    """Generates a unique identifier string for a guide using its genomic and
    sequence attributes.

    The identifier includes chromosome, start, stop, strand, sgRNA, and PAM sequence.

    Args:
        chrom: Chromosome name.
        start: Start position of the guide.
        stop: Stop position of the guide.
        strand: Strand orientation ('+' or '-').
        sgrna: The sgRNA sequence.
        pam: The PAM sequence.

    Returns:
        str: A formatted string uniquely identifying the guide.
    """
    return f"{chrom}_{start}_{stop}_{strand}_{sgrna}_{pam}"


def _assign_extended_guide_ids(report: pd.DataFrame) -> pd.DataFrame:
    """Adds a 'guide_id' column to the report DataFrame using chromosome, start,
    and strand.

    This function generates a unique identifier for each guide and appends it as
    a new column.

    Args:
        report: DataFrame containing guide information.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'guide_id' column.
    """
    # compute guide ids and drop non unique sites
    report["guide_id"] = report.apply(
        lambda x: _compute_extended_guide_id(
            x[REPORTCOLS[0]],
            x[REPORTCOLS[1]],
            x[REPORTCOLS[2]],
            x[REPORTCOLS[6]],
            x[REPORTCOLS[3]],
            x[REPORTCOLS[4]],
        ),
        axis=1,
    )
    report = report.drop_duplicates(subset="guide_id")
    return report


def _compute_guide_id(chrom: str, start: int, strand: str) -> str:
    """Creates a unique identifier string for a guide using chromosome, start
    position, and strand.

    The identifier is useful for referencing guides in downstream analyses and
    grouping operations.

    Args:
        chrom: Chromosome name.
        start: Start position of the guide.
        strand: Strand orientation ('+' or '-').

    Returns:
        str: A formatted string uniquely identifying the guide.
    """
    return f"{chrom}_{start}_{strand}"


def _assign_guide_ids(report: pd.DataFrame) -> pd.DataFrame:
    """Adds a 'guide_id' column to the report DataFrame using chromosome, start
    position, and strand.

    This function generates a unique identifier for each guide and appends it as
    a new column.

    Args:
        report: DataFrame containing guide information.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'guide_id' column.
    """
    report["guide_id"] = report.apply(
        lambda x: _compute_guide_id(
            x[REPORTCOLS[0]], x[REPORTCOLS[1]], x[REPORTCOLS[6]]
        ),
        axis=1,
    )
    return report


def _assess_guide_type(origin: str, sgrna: str, pam: str, debug: bool) -> int:
    """Assesses a guide type based on its origin, sgRNA, and PAM sequence.

    Returns an integer representing the guide type: reference, spacer+PAM alternative,
    spacer alternative, or PAM alternative.

    Args:
        origin: The origin of the guide ('ref' for reference, otherwise alternative).
        sgrna: The sgRNA sequence.
        pam: The PAM sequence.
        debug: Boolean flag for debug mode.

    Returns:
        int: Guide type (0 for reference, 1 for spacer+PAM alternative, 2 for
            spacer alternative, 3 for PAM alternative).

    Raises:
        CrisprHawkGraphicalReportsError: If the guide type cannot be determined.
    """
    # types: 0 -> ref; 1 -> spacer+pam alt; 2 -> spacer alt; 3 -> pam alt
    if origin == "ref":  # reference type guide
        return 0
    assert origin != "ref"  # asses alternative guide type
    sgrna_alt = is_lowercase(sgrna)  # check sgrna
    pam_alt = is_lowercase(pam)  # check pam
    if sgrna_alt:  # spacer alt or spacer+pam alt
        return 1 if pam_alt else 2
    if pam_alt:  # pam alt
        return 3
    # this chunk of code should never be reached
    exception_handler(
        CrisprHawkGraphicalReportsError,
        f"Unknown guide type for grna {sgrna}",
        os.EX_DATAERR,
        debug,
    )


def _assign_guide_type(report: pd.DataFrame, debug: bool) -> pd.DataFrame:
    """Adds a 'guide_id' column to the report DataFrame using chromosome, start,
    and strand.

    This function generates a unique identifier for each guide and appends it as
    a new column.

    Args:
        report: DataFrame containing guide information.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'guide_id' column.
    """
    report["guide_type"] = report.apply(
        lambda x: _assess_guide_type(x["origin"], x["sgRNA_sequence"], x["pam"], debug),
        axis=1,
    )  # assign guide type
    return report


def _count_guide_type(guide_types: List[int]) -> Dict[str, int]:
    """Counts the number of guides for each guide type and returns a summary
    dictionary.

    The function returns a dictionary mapping guide type labels to their respective
    counts.

    Args:
        guide_types: A list of integers representing guide types.

    Returns:
        Dict[str, int]: A dictionary with guide type labels as keys and counts
            as values.
    """
    types_data = {label: 0 for _, label in GUIDETYPES.items()}
    for gt in guide_types:  # count number of guides for each type
        types_data[GUIDETYPES[gt]] += 1
    return types_data


def _draw_piechart(
    data: Dict[str, int], region_format: str, prefix: str, outdir: str
) -> None:
    """Draws and saves a pie chart visualizing the distribution of guide types
    for a region.

    The pie chart is saved as a PNG file in the specified output directory.

    Args:
        data: Dictionary mapping guide type labels to their counts.
        region_format: String representing the region for labeling the chart.
        prefix: Prefix for the output file name.
        outdir: Directory where the pie chart image will be saved.

    Returns:
        None
    """
    labels = list(data.keys())  # pie chart labels
    values = list(data.values())  # pie chart data
    colors = ["#5f8dd3ff", "#0055d4ff", "#ff6600ff", "#ffcc00ff"]
    explode = (0, 0, 0.05, 0.1)
    f, ax = plt.subplots(1, 1, figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(
        values,
        explode=explode,
        colors=colors,
        autopct="%.2f%%",
        shadow=False,
        startangle=140,
        textprops={"fontsize": 14},
        pctdistance=1.1,
    )
    ax.set_title(f"Guide Types (region: {region_format})", fontsize=20)
    plt.legend(labels, loc=(0.8, 0), prop={"size": 16})
    plt.axis("equal")
    plt.tight_layout()
    piechart_fname = os.path.join(outdir, f"{prefix}_guides_type.png")
    plt.savefig(piechart_fname, format="png", dpi=300)


def piechart_guides_type(
    report: pd.DataFrame,
    region: Region,
    prefix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Generates and saves a pie chart showing the distribution of guide types
    for a genomic region.

    The function processes the report DataFrame, assigns guide types, and visualizes
    their proportions in a pie chart.

    Args:
        report: DataFrame containing guide information.
        region: Region object representing the genomic region.
        prefix: Prefix for the output file name.
        outdir: Directory where the pie chart image will be saved.
        verbosity: Verbosity level for logging.
        debug: Boolean flag for debug mode.

    Returns:
        None
    """
    print_verbosity("Computing guide type pie chart", verbosity, VERBOSITYLVL[3])
    start = time()
    report = _assign_extended_guide_ids(report)  # compute guide ids
    report = _assign_guide_type(report, debug)  # assign guide type
    # draw pie chart for guide types
    _draw_piechart(
        _count_guide_type(report["guide_type"].tolist()),
        str(region.coordinates),
        prefix,
        outdir,
    )
    print_verbosity(
        f"Guide type pie chart computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )


def _compute_nsamples(samples: str) -> int:
    """Computes the number of unique samples from a sample string.

    Returns 0 if the value is 'REF', otherwise returns the count of unique sample
    names.

    Args:
        samples: A string representing sample information, with samples separated
            by commas.

    Returns:
        int: The number of unique samples.
    """
    if samples == "REF":
        return 0
    return len({s.split(":")[0] for s in samples.split(",")})


def _assign_nsamples(report: pd.DataFrame) -> pd.DataFrame:
    """Adds a 'n_samples' column to the report DataFrame by computing the number
    of unique samples for each row.

    This function calculates the sample count for each guide and appends it as
    a new column.

    Args:
        report: DataFrame containing guide information.

    Returns:
        pd.DataFrame: The input DataFrame with an added 'n_samples' column.
    """
    report["n_samples"] = report.apply(
        lambda x: _compute_nsamples(x[REPORTCOLS[14]]), axis=1
    )
    return report


def _check_candidate_guides(
    region: Region, candidate_guides: List[CandidateGuide], report: pd.DataFrame
) -> List[str]:
    """Filters and returns candidate guide IDs that are contained within a given
    region.

    This function checks each candidate guide and includes it if its coordinate
    is within the region's coordinates.


    Args:
        region: Region object representing the genomic region of interest.
        candidate_guides: List of CandidateGuide objects to be checked.
        report: DataFrame containing guide information (not used in this function).

    Returns:
        List[str]: List of guide ID strings for candidate guides within the region.
    """
    return [
        f"{cg.contig}_{cg.position}_{cg.strand}"
        for cg in candidate_guides
        if region.coordinates.contains(cg.coordinate)
    ]


def _compute_group_delta(group: Any, score: str) -> pd.Series:
    """Calculates the delta score for each guide in a group relative to the
    reference guide.

    This function subtracts the reference guide's score from each guide's score
    within the group.

    Args:
        group: DataFrame containing guides for a single guide ID.
        score: The score column to use for delta calculation.

    Returns:
        pd.Series: The group DataFrame with an updated delta column.
    """
    refrow = group[group[REPORTCOLS[13]] == "ref"]  # reference grna
    if refrow.empty:  # only alternative grna
        return group
    refscore = refrow[score].values[0]
    group[DELTACOLS[0]] = group[score] - refscore  # compute delta between scores
    return group


def _compute_deltas(report: pd.DataFrame, score: str) -> pd.DataFrame:
    """Calculates delta values for each guide in the report relative to the reference
    guide.

    This function initializes the delta column and applies group-wise delta
    computation for each guide ID.

    Args:
        report: DataFrame containing guide information.
        score: The score column to use for delta calculation.

    Returns:
        pd.DataFrame: The input report DataFrame with an updated delta column.
    """
    report[DELTACOLS[0]] = 0.0  # initialize deltas
    return report.groupby("guide_id", group_keys=False).apply(
        _compute_group_delta, score
    )


def _compute_group_delta_abs(group: Any, score: str) -> pd.Series:
    """Calculates both the delta and absolute delta scores for each guide in a
    group relative to the reference guide.

    This function computes the difference and absolute difference between each
    guide's score and the reference guide's score within the group.

    Args:
        group: DataFrame containing guides for a single guide ID.
        score: The score column to use for delta calculation.

    Returns:
        pd.Series: The group DataFrame with updated delta and absolute delta columns.
    """
    refrow = group[group[REPORTCOLS[13]] == "ref"]  # reference grna
    if refrow.empty:  # only alternative grna
        return group
    refscore = refrow[score].values[0]
    group[DELTACOLS[0]] = group[score] - refscore
    group[DELTACOLS[1]] = np.abs(group[DELTACOLS[0]])  # compute abs delta
    return group


def _compute_deltas_abs(report: pd.DataFrame, score: str) -> pd.DataFrame:
    """Calculates both delta and absolute delta values for each guide in the
    report relative to the reference guide.

    This function initializes the delta and absolute delta columns, then applies
    group-wise computation for each guide ID.

    Args:
        report: DataFrame containing guide information.
        score: The score column to use for delta calculation.

    Returns:
        pd.DataFrame: The input report DataFrame with updated delta and absolute
            delta columns.
    """
    # initialize deltas and absolute deltas
    report[DELTACOLS[0]] = 0.0
    report[DELTACOLS[1]] = 0.0
    return report.groupby("guide_id", group_keys=False).apply(
        _compute_group_delta_abs, score
    )


def _compute_scores_delta(report: pd.DataFrame, score: str) -> pd.DataFrame:
    """Computes delta or absolute delta values for each guide in the report based
    on the score type.

    This function selects the appropriate delta computation method depending on
    the score column.

    Args:
        report: DataFrame containing guide information.
        score: The score column to use for delta calculation.

    Returns:
        pd.DataFrame: The input report DataFrame with updated delta or absolute
            delta columns.
    """
    if score in SCORES[3:]:  # cfdon, elevationon
        return _compute_deltas(report, score)
    return _compute_deltas_abs(report, score)  # azimuth, rs3, deepcpf1


def _filter_valid_alts(alts: pd.DataFrame, refscore: float, score: str):
    """Filters alternative guides based on their score relative to the reference
    score.

    For certain score types, only alternatives with a lower score than the reference
    are included; otherwise, all alternatives are returned.

    Args:
        alts: DataFrame containing alternative guide information.
        refscore: The score of the reference guide.
        score: The score column to use for filtering alternatives.

    Returns:
        pd.DataFrame: DataFrame of valid alternative guides.
    """
    return alts[alts[score] < refscore] if score in SCORES[3:] else alts.copy()


def _extract_alt_data(alt_row: pd.Series, score: str) -> Dict[str, Any]:
    """Extracts alternative guide information from a DataFrame row as a dictionary.

    The function returns a dictionary containing sgRNA, PAM, score, delta, absolute
    delta, sample count, and variant ID for an alternative guide.

    Args:
        alt_row: A pandas Series representing a row with alternative guide data.
        score: The score column to extract for the alternative guide.

    Returns:
        Dict[str, Any]: Dictionary with alternative guide information.
    """
    return {
        "alt_sgRNA": alt_row["sgRNA_sequence"],
        "pam": alt_row["pam"],
        "alt_score": alt_row[score],
        "delta": alt_row["delta"],
        "abs_delta": alt_row.get("abs_delta", abs(alt_row["delta"])),
        "n_samples": alt_row["n_samples"],
        "variant_id": alt_row["variant_id"],
    }


def _build_guide_rows(report: pd.DataFrame, score: str) -> Dict[str, Any]:
    """Builds a dictionary of guide data for each guide ID, including reference
    and alternative guides.

    The function groups the DataFrame by guide ID and processes each group to
    extract relevant guide information.

    Args:
        report: DataFrame containing guide information.
        score: The score column to use for processing guides.

    Returns:
        Dict[str, Any]: A dictionary mapping guide IDs to their reference and
            alternative guide data.
    """
    grouped = report.groupby("guide_id", sort=False)
    guide_rows = {}
    for guide_id, group in grouped:
        ref = group[group[REPORTCOLS[13]] == "ref"]
        alts = group[group[REPORTCOLS[13]] == "alt"]
        if ref.empty:
            continue
        refscore = ref[score].values[0]
        refseq = ref[REPORTCOLS[3]].values[0]
        refpam = ref[REPORTCOLS[4]].values[0]
        refsamples = ref["n_samples"].values[0]
        # filter and extract alternative data
        valid_alts = _filter_valid_alts(alts, refscore, score)
        alt_list = [
            _extract_alt_data(alt_row, score) for _, alt_row in valid_alts.iterrows()
        ]
        guide_rows[guide_id] = {
            "ref_sgRNA": refseq,
            "pam": refpam,
            "ref_score": refscore,
            "ref_samples": refsamples,
            "alts": alt_list,
        }
    return guide_rows


def _calculate_worst_delta(alts: Dict[str, Any], score: str) -> float:
    """Calculates the worst delta value among alternative guides for a given
    score type.

    For certain score types, returns the maximum absolute delta; for others,
    returns the minimum delta.

    Args:
        alts: List of dictionaries containing alternative guide information.
        score: The score column to use for determining the worst delta.

    Returns:
        float: The worst delta value among the alternatives.
    """
    if not alts:
        return 0.0
    if score in SCORES[:3]:
        return max(alt[DELTACOLS[1]] for alt in alts)  # type: ignore
    return min(alt[DELTACOLS[0]] for alt in alts)  # type: ignore


def _rank_guides_by_worst_delta(guide_rows: Dict[str, Any], score: str) -> pd.DataFrame:
    """Ranks guides by their worst delta value and returns a sorted DataFrame.

    This function computes the worst delta for each guide and sorts the guides
    in ascending or descending order depending on the score type.

    Args:
        guide_rows: Dictionary mapping guide IDs to their reference and alternative
            guide data.
        score: The score column to use for ranking guides.

    Returns:
        pd.DataFrame: DataFrame of guides sorted by their worst delta value.
    """
    worst_alts = []
    for gid, data in guide_rows.items():
        worst_delta = _calculate_worst_delta(data["alts"], score)
        worst_alts.append((gid, worst_delta))
    return (
        pd.DataFrame(worst_alts, columns=["guide_id", "delta"])
        .sort_values("delta", ascending=(score in SCORES[3:]))
        .reset_index(drop=True)
    )


def _mask(df: pd.DataFrame, mask: str) -> pd.Series:
    """Creates a boolean mask for selecting rows in a DataFrame by guide ID.

    This function returns a boolean Series indicating which rows have a guide_id
    matching the provided mask.

    Args:
        df: DataFrame containing a 'guide_id' column.
        mask: The guide ID to match.

    Returns:
        pd.Series: Boolean Series where True indicates a matching guide ID.
    """
    return df["guide_id"] == mask


def _candidate_guides_subreports(
    worst_df: pd.DataFrame, cgids: List[str]
) -> pd.DataFrame:
    """Selects candidate guides from the ranked DataFrame and fills remaining
    slots with top non-candidate guides.

    This function returns a DataFrame containing all specified candidate guides
    and the highest-ranked non-candidate guides up to a total of 25.

    Args:
        worst_df: DataFrame containing guides ranked by their worst delta values.
        cgids: List of candidate guide IDs to include in the selection.

    Returns:
        pd.DataFrame: DataFrame containing candidate guides and additional top guides.
    """
    cgrows = pd.concat(
        [pd.DataFrame([worst_df[_mask(worst_df, cgid)].iloc[0]]) for cgid in cgids]
    )
    cgmask_ = ~worst_df["guide_id"].isin(cgids)
    others = worst_df[cgmask_].head(25 - len(cgids))
    return pd.concat([cgrows, others])


def _select_top_guides(worst_df: pd.DataFrame, cgids: List[str]) -> pd.DataFrame:
    """Selects the top guides for reporting, prioritizing candidate guides if provided.

    This function returns a DataFrame of the top 25 guides, including all candidate
    guides if specified, and assigns a rank to each.

    Args:
        worst_df: DataFrame containing guides ranked by their worst delta values.
        cgids: List of candidate guide IDs to prioritize in the selection.

    Returns:
        pd.DataFrame: DataFrame of the top 25 guides with assigned ranks.
    """
    if cgids:  # canidate guide ids
        guidesfinal = _candidate_guides_subreports(worst_df, cgids)
    else:
        guidesfinal = worst_df.head(25)
    guidesfinal = guidesfinal.reset_index(drop=True)
    guidesfinal["Rank"] = guidesfinal.index + 1  # add rank column (1-indexed)
    return guidesfinal


def _fill_delta_table_row_data(
    outrow: Dict[str, Any], alt: pd.Series, idx: int
) -> Dict[str, Any]:
    """Fills a delta table row dictionary with alternative guide data for a
    specific index.

    This function updates the output row with alternative guide information such
    as sgRNA, PAM, score, delta, absolute delta, sample count, and variant ID.

    Args:
        outrow: Dictionary representing a row in the delta table.
        alt: Series or dictionary containing alternative guide data.
        idx: The index of the alternative guide (1-based).

    Returns:
        Dict[str, Any]: The updated output row dictionary with alternative guide data.
    """
    outrow[f"alt{idx}_sgRNA"] = alt["alt_sgRNA"]
    outrow[f"alt{idx}_pam"] = alt["pam"]
    outrow[f"alt{idx}_score"] = alt["alt_score"]
    outrow[f"alt{idx}_delta"] = alt["delta"]
    outrow[f"alt{idx}_abs_delta"] = alt["abs_delta"]
    outrow[f"alt{idx}_n_samples"] = alt["n_samples"]
    outrow[f"alt{idx}_variant_id"] = alt["variant_id"]
    return outrow


def _fill_delta_table_row_nan(outrow: Dict[str, Any], idx: int) -> Dict[str, Any]:
    """Fills a delta table row dictionary with NaN values for a missing alternative
    guide.

    This function updates the output row with NaN for all alternative guide fields
    at the specified index.

    Args:
        outrow: Dictionary representing a row in the delta table.
        idx: The index of the alternative guide (1-based).

    Returns:
        Dict[str, Any]: The updated output row dictionary with NaN values for the
            alternative guide.
    """
    outrow[f"alt{idx}_sgRNA"] = np.nan
    outrow[f"alt{idx}_pam"] = np.nan
    outrow[f"alt{idx}_score"] = np.nan
    outrow[f"alt{idx}_delta"] = np.nan
    outrow[f"alt{idx}_abs_delta"] = np.nan
    outrow[f"alt{idx}_n_samples"] = np.nan
    outrow[f"alt{idx}_variant_id"] = np.nan
    return outrow


def _build_delta_table_row(
    guide_id: str, rank: int, data: Any, max_alts: int
) -> Dict[str, Any]:
    """Builds a dictionary representing a row in the delta table for a guide and
    its alternatives.

    The function creates a row with reference guide information and adds columns
    for each alternative guide, filling with NaN if alternatives are missing.

    Args:
        guide_id: The unique identifier for the guide.
        rank: The rank of the guide.
        data: Dictionary containing reference and alternative guide data.
        max_alts: Maximum number of alternative guides to include.

    Returns:
        Dict[str, Any]: Dictionary representing a row in the delta table with
            reference and alternative guide information.
    """
    outrow = {
        "guide_id": guide_id,
        "Rank": rank,
        "ref_sgRNA": data["ref_sgRNA"],
        "pam": data["pam"],
        "ref_score": data["ref_score"],
        "ref_n_samples": data["ref_samples"],
    }
    for i in range(max_alts):  # add columns for each alternative (nan if not present)
        alt_num = i + 1
        if i < len(data["alts"]):
            alt = data["alts"][i]
            outrow = _fill_delta_table_row_data(outrow, alt, alt_num)
        else:  # fill with NaN for missing alternatives
            outrow = _fill_delta_table_row_nan(outrow, alt_num)
    return outrow


def _construct_delta_table(
    guidesfinal: pd.DataFrame, guide_rows: Dict[str, Any]
) -> pd.DataFrame:
    """Constructs a wide-format DataFrame summarizing reference and alternative
    guide information for top-ranked guides.

    The function returns a DataFrame with columns for guide ID, rank, reference
    guide details, and alternative guide details.

    Args:
        guidesfinal: DataFrame containing the top-ranked guides and their ranks.
        guide_rows: Dictionary mapping guide IDs to their reference and alternative
            guide data.

    Returns:
        pd.DataFrame: Wide-format DataFrame summarizing guide and alternative information.
    """
    # determine maximum number of alternatives across all selected guides
    max_alts = max(len(guide_rows[gid]["alts"]) for gid in guidesfinal["guide_id"])
    rows = []
    for _, row in guidesfinal.iterrows():
        gid = row["guide_id"]
        data = guide_rows[gid]
        out_row = _build_delta_table_row(gid, row["Rank"], data, max_alts)
        rows.append(out_row)
    return pd.DataFrame(rows)


def _compute_delta_table(
    report: pd.DataFrame, cgids: List[str], score: str
) -> pd.DataFrame:
    """Computes a summary table of delta scores for guides and their alternatives.

    This function processes the report DataFrame, calculates delta values, ranks guides,
    selects the top guides, and constructs a wide-format table summarizing reference and
    alternative guide information.

    Args:
        report: DataFrame containing guide and score information.
        cgids: List of candidate guide IDs to prioritize in the selection.
        score: The score column to use for delta calculation.

    Returns:
        pd.DataFrame: Wide-format DataFrame summarizing delta scores for top
            guides and their alternatives.
    """
    report = _assign_guide_ids(report)  # add guide ids for grouping
    report = _assign_nsamples(report)  # compute number of samples
    report_deltas = _compute_scores_delta(report, score)  # compute deltas
    guide_rows = _build_guide_rows(report_deltas, score)  # build guide data structure
    worst_df = _rank_guides_by_worst_delta(guide_rows, score)  # rank by worst delta
    guidesfinal = _select_top_guides(worst_df, cgids)  # select top guides
    return _construct_delta_table(guidesfinal, guide_rows)  # construct output table


def _save_plot(output_prefix: str, score_col: str, output_dir: str) -> str:
    """Saves the current plot to a PNG file in the specified output directory.

    The function returns the path to the saved plot file.

    Args:
        output_prefix: Prefix for the output file name.
        score_col: The score column used in the plot, included in the file name.
        output_dir: Directory where the plot image will be saved.

    Returns:
        str: Path to the saved plot file.
    """
    filename = f"{output_prefix}_{score_col}_delta.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, format="png", dpi=DPI, bbox_inches="tight")
    return output_path


def _compute_rank_id(rank: int, guideid: str) -> str:
    """Creates a combined string of rank and genomic coordinate for a guide.

    This function returns a formatted string with the guide's rank and its
    chromosome and position.

    Args:
        rank: The rank of the guide.
        guideid: The unique identifier for the guide, containing chromosome and
            position.

    Returns:
        str: A formatted string with rank and genomic coordinate.
    """
    chrom, pos = guideid.split("_")[:2]
    return f"Rank {rank}, {chrom}:{pos}"


def _prepare_top_guides(df: pd.DataFrame) -> pd.DataFrame:
    """Prepares a DataFrame of the top 25 guides sorted by rank and adds a combined
    rank and genomic coordinate column.

    This function sorts the DataFrame by the 'Rank' column, selects the top 25 guides,
    and adds a 'rank_chr_start' column for plotting.

    Args:
        df: DataFrame containing guide information with 'Rank' and 'guide_id'
            columns.

    Returns:
        pd.DataFrame: DataFrame of the top 25 guides with an added 'rank_chr_start'
            column.
    """
    colnames = df.columns.tolist()
    assert "Rank" in colnames and "guide_id" in colnames
    top_df = df.sort_values("Rank").head(25).copy()
    top_df["rank_chr_start"] = top_df.apply(
        lambda x: _compute_rank_id(x["Rank"], x["guide_id"]), axis=1
    )
    return top_df


def _has_candidate_guide(df: pd.DataFrame, cgids: List[str]) -> List[str]:
    """Identifies candidate guide IDs that are present in the DataFrame.

    This function returns a list of candidate guide IDs that exist in the 'guide_id'
    column of the DataFrame.

    Args:
        df: DataFrame containing a 'guide_id' column.
        cgids: List of candidate guide IDs to check for presence.

    Returns:
        List[str]: List of candidate guide IDs found in the DataFrame.
    """
    guideids = set(df["guide_id"].tolist())
    return [cgid for cgid in cgids if cgid if cgid in guideids]


def _has_variant(row: pd.Series, colnames: List[str]) -> bool:
    """Determines if a row contains at least one alternative variant in the
    specified columns.

    This function checks if any of the given columns in the row are not NaN and
    not equal to 'REF'.

    Args:
        row: A pandas Series representing a row of guide and variant data.
        colnames: List of column names to check for alternative variants.

    Returns:
        bool: True if at least one alternative variant is present, False otherwise.
    """
    return any(pd.notna(row[c]) and row[c] != "REF" for c in colnames)


def _extract_variant_keys(df: pd.DataFrame) -> List[str]:
    """Extracts unique variant keys from a DataFrame for guides with alternative
    variants.

    This function returns a list of 'rank_chr_start' values for rows that have
    at least one alternative variant.

    Args:
        df: DataFrame containing guide and alternative variant information.

    Returns:
        List[str]: List of unique variant keys for guides with alternatives.
    """
    colnames = df.columns.tolist()
    alt_samples_col = [
        col for col in colnames if col.startswith("alt") and col.endswith("_n_samples")
    ]
    variant_keys = {
        row["rank_chr_start"]
        for _, row in df.iterrows()
        if _has_variant(row, alt_samples_col)
    }
    return list(variant_keys)


def _generate_color_palette(n_colors: int) -> List[Tuple[float, float, float]]:
    """Generates a color palette with a specified number of distinct colors for
    plotting.

    This function creates a list of RGB color tuples by sampling from predefined
    color maps, ensuring enough unique colors for visualization.

    Args:
        n_colors: The number of distinct colors to generate.

    Returns:
        List[Tuple[float, float, float]]: List of RGB color tuples.
    """
    # Start with darkest shade (index 7) from each colormap
    palette = [sns.color_palette(cmap, 9)[7] for cmap in BASE_CMAPS]
    # if we need more colors, add medium-dark shades (indices 6-8)
    if n_colors > len(palette):
        extra_colors = []
        for shade_idx in range(6, 9):
            extra_colors.extend(
                sns.color_palette(cmap, 9)[shade_idx] for cmap in BASE_CMAPS
            )
        needed = n_colors - len(palette)
        palette += extra_colors[:needed]
    random.shuffle(palette)
    return palette


def _assign_variant_colors(
    variant_keys: List[str],
) -> Dict[str, Tuple[float, float, float]]:
    """Assigns a unique color to each variant key for plotting purposes.

    This function generates a color palette and maps each variant key to a color.

    Args:
        variant_keys: List of unique variant keys to assign colors to.

    Returns:
        Dict[str, Tuple[float, float, float]]: Dictionary mapping variant keys to
            RGB color tuples.
    """
    palette = _generate_color_palette(len(variant_keys))
    return {key: palette[i] for i, key in enumerate(variant_keys)}


def _get_legend_size(nsamples: int, debug: bool) -> int:
    """Determines the marker size for plotting based on the number of samples.

    This function returns a scaled marker size for use in legends and scatter
    plots, with larger sizes for higher sample counts. If the number of samples
    is invalid, an exception is raised.

    Args:
        nsamples: The number of samples represented by a guide.
        debug: Boolean flag for debug mode, used in exception handling.

    Returns:
        int: The calculated marker size for plotting.

    Raises:
        CrisprHawkGraphicalReportsError: If the number of samples is invalid.
    """
    thresholds = [1, 20, 50, 100, 200, np.inf]
    bases = [1, 10, 35, 75, 150, 300]
    for limit, base in zip(thresholds, bases):
        if nsamples <= limit:
            return 150 * np.sqrt(base)
    exception_handler(
        CrisprHawkGraphicalReportsError,
        f"Invalid number of samples ({nsamples})",
        os.EX_DATAERR,
        debug,
    )


def _get_marker_style(nsamples: int) -> str:
    """Determines the marker style for plotting based on the number of samples.

    This function returns a diamond marker for single-sample guides and a circle
    marker for guides with more than one sample.

    Args:
        nsamples: The number of samples represented by a guide.

    Returns:
        str: The marker style ('D' for diamond, 'o' for circle).
    """
    return "D" if nsamples == 1 else "o"


def _plot_alternative_guides(
    ax: Axes,
    df: pd.DataFrame,
    palette: Dict[str, Tuple[float, float, float]],
    debug: bool,
) -> None:
    """Plots alternative guide scores as scatter points on the provided axes.

    This function iterates through the DataFrame rows and plots each alternative
    guide's score, using color and marker size to represent variant and sample
    count, respectively.

    Args:
        ax: Matplotlib Axes object to plot on.
        df: DataFrame containing guide and alternative score information.
        palette: Dictionary mapping variant keys to RGB color tuples.

    Returns:
        None
    """
    for _, row in df.iterrows():
        rank, rank_id = row["Rank"], row["rank_chr_start"]
        color = palette.get(rank_id, "gray")
        for i in range(1, 1000):  # iterate through alternative columns
            score_cname = f"alt{i}_score"
            nsamples_cname = f"alt{i}_n_samples"
            if any(c not in row for c in [score_cname, nsamples_cname]):
                break  # stop when run out of alternative guides columns
            if pd.isna(row[score_cname]):  # skip if no data for thi alt guide
                continue
            # get sample count and corresponding marker properties
            nsamples = int(row[nsamples_cname]) if pd.notna(row[nsamples_cname]) else 1
            size = _get_legend_size(nsamples, debug)
            marker = _get_marker_style(nsamples)
            ax.scatter(rank, row[score_cname], color=color, s=size, alpha=0.6, edgecolors="white", linewidth=0.5, marker=marker, zorder=3)  # type: ignore # plot guide


def _plot_reference_line(ax: Axes) -> Line2D:
    """Plots a horizontal reference line at y=1 on the provided axes.

    This function adds a gray reference line to the plot and returns a legend
    handle for the reference.

    Args:
        ax: Matplotlib Axes object to plot on.

    Returns:
        Line2D: The legend handle for the reference line.
    """
    ax.axhline(y=1, color="gray", linestyle="-", alpha=0.6, linewidth=3, zorder=5)
    return Line2D([0], [0], color="gray", linewidth=3, label="Reference", alpha=0.6)


def _plot_reference_points(ax: Axes, df: pd.DataFrame) -> Line2D:
    """Plots reference guide scores as points on the provided axes.

    This function adds black scatter points for reference guides and returns a
    legend handle for the reference.

    Args:
        ax: Matplotlib Axes object to plot on.
        df: DataFrame containing guide information with 'Rank' and 'ref_score' columns.

    Returns:
        Line2D: The legend handle for the reference points.
    """
    ax.scatter(
        df["Rank"],
        df["ref_score"],
        color="black",
        s=600,
        zorder=5,
        edgecolors="white",
        linewidth=0.7,
        label="Reference",
        alpha=0.8,
    )
    return Line2D(
        [0],
        [0],
        marker="o",
        color="w",
        markerfacecolor="dimgray",
        markersize=np.sqrt(600),
        markeredgecolor="white",
        label="Reference",
        alpha=0.8,
    )


def _plot_reference_elements(ax: Axes, df: pd.DataFrame, score: str) -> Line2D:
    """Plots reference guide elements (line or points) on the provided axes based
    on the score type.

    This function adds either a horizontal reference line or reference points to
    the plot, depending on the score column.

    Args:
        ax: Matplotlib Axes object to plot on.
        df: DataFrame containing guide information.
        score: The score column used to determine the reference plot style.

    Returns:
        Line2D: The legend handle for the reference element.
    """
    if score in SCORES[3:]:  # cfdon and elevationon
        return _plot_reference_line(ax)  # reference guides as line on top
    return _plot_reference_points(ax, df)  # reference guides as dots


def _configure_x_axis(ax: Axes, df: pd.DataFrame, cgdids_present: List[str]) -> None:
    """Configures the x-axis of the plot with guide labels and highlights
    candidate guides.

    This function sets the x-tick labels to guide sequences, rotates them for
    readability, and bolds the labels for candidate guides if provided.

    Args:
        ax: Matplotlib Axes object to configure.
        df: DataFrame containing guide information with 'Rank' and 'ref_sgRNA' columns.
        cgdids_present: List of candidate guide IDs to highlight in the labels.

    Returns:
        None
    """
    xtick_labels = [row["ref_sgRNA"] for _, row in df.iterrows()]
    ax.set_xticks(df["Rank"])
    ax.set_xticklabels(xtick_labels, rotation=45, ha="right", fontsize=15)
    if cgdids_present:  # candidate guides requested
        for label, (guideid, _) in zip(
            ax.get_xticklabels(), df[["guide_id", "Rank"]].values
        ):
            if guideid in cgdids_present:
                label.set_fontweight("bold")


def _get_axis_labels_and_title(score: str) -> Tuple[str, str, str]:
    """Generates axis labels and plot title based on the score type.

    This function returns the x-axis label, y-axis label, and plot title appropriate
    for the given score column.

    Args:
        score: The score column used to determine axis labels and title.

    Returns:
        Tuple[str, str, str]: The x-axis label, y-axis label, and plot title.
    """
    xlabel = "Guide"
    if score in SCORES[3:]:  # cfdon and elevationon
        ylabel = "Variant Effect (CFD)"
        title = "Variant Effect on Alternative On-Targets"
    else:
        score_name = score.split("_")[1].upper()
        ylabel = f"On-Target Efficiency ({score_name})"
        title = "Guides On-Target Efficiency"
    return xlabel, ylabel, title


def _configure_y_axis(ax: Axes, score: str) -> None:
    """Configures the y-axis of the plot with appropriate limits, labels, and
    title based on the score type.

    This function sets the y-axis range, axis labels, and plot title for the
    given score column.

    Args:
        ax: Matplotlib Axes object to configure.
        score: The score column used to determine y-axis limits and labels.

    Returns:
        None
    """
    ax.tick_params(axis="y", labelsize=15)
    if score == SCORES[1]:  # rs3
        ax.set_ylim(-2.05, 2.05)
    elif score == SCORES[2]:  # deepcpf1
        ax.set_ylim(-10, 110)
    else:
        ax.set_ylim(-0.05, 1.05)
    # get and set labels
    xlabel, ylabel, title = _get_axis_labels_and_title(score)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title, fontsize=24)


def _create_size_legend_handles(ref_legend_handle: Line2D) -> List[Line2D]:
    """Creates legend handles for different sample sizes and includes the
    reference guide handle.

    This function generates a list of Line2D objects for the legend, representing
    marker sizes for various sample counts and the reference guide.

    Args:
        ref_legend_handle: The legend handle for the reference guide.

    Returns:
        List[Line2D]: List of legend handles for sample sizes and the reference guide.
    """
    scaled_sizes = [150 * np.sqrt(n) for n in LEGEND_SAMPLE_COUNTS]
    size_handles = [
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
        for label, size, marker in zip(LEGEND_LABELS, scaled_sizes, LEGEND_MARKERS)
    ]
    size_handles.insert(0, ref_legend_handle)  # add reference legend at the beginning
    return size_handles


def _add_legend(ref_legend_handle: Line2D) -> None:
    """Adds a legend to the plot showing the number of samples for each guide marker.

    This function creates and configures a legend with marker size and style
    information, including a title.

    Args:
        ref_legend_handle: The legend handle for the reference guide.

    Returns:
        None
    """
    handles = _create_size_legend_handles(ref_legend_handle)
    legend = plt.legend(
        handles=handles,
        title="Number of samples",
        frameon=False,
        bbox_to_anchor=(0.5, -0.55),
        loc="lower center",
        ncol=7,
        fontsize=11,
        handletextpad=1,
        columnspacing=5,
        labelspacing=2,
    )
    legend.get_title().set_fontsize(24)


def _apply_plot_styling(ax: Axes) -> None:
    """Applies grid and styling to the plot axes for improved visualization.

    This function adds a light dashed grid and removes the top and right spines
    for a cleaner plot appearance.

    Args:
        ax: Matplotlib Axes object to style.

    Returns:
        None
    """
    ax.grid(True, alpha=0.3, linestyle="--")
    sns.despine()


def _adjust_layout(f: Figure) -> None:
    """Adjusts the layout of the figure to provide additional space at the bottom.

    This function modifies the subplot parameters to prevent label overlap and
    improve plot readability.

    Args:
        f: Matplotlib Figure object to adjust.

    Returns:
        None
    """
    f.subplots_adjust(bottom=0.25)


def _dotplot_delta(
    df: pd.DataFrame,
    score: str,
    cgids: List[str],
    prefix: str,
    outdir: str,
    debug: bool,
) -> None:
    """Creates and saves a delta dot plot visualizing guide efficiency scores and
    their alternatives.

    This function prepares the data, generates a dot plot for reference and alternative
    guides, applies styling and legends, and saves the plot to the specified output
    directory.

    Args:
        df: DataFrame containing guide and alternative score information.
        score: The score column to visualize.
        cgids: List of candidate guide IDs to highlight in the plot.
        prefix: Prefix for the output file name.
        outdir: Directory where the plot image will be saved.

    Returns:
        None
    """
    topdf = _prepare_top_guides(df)  # prepare data for delta plots
    cgids_present = _has_candidate_guide(
        df, cgids
    )  # check if candidate guides provided
    variant_keys = _extract_variant_keys(topdf)  # retrieve alt guides rank ids
    variant_colors = _assign_variant_colors(variant_keys)  # generate color palette
    f, ax = plt.subplots(1, 1, figsize=FIGURE_SIZE)  # initialize figure
    _plot_alternative_guides(ax, topdf, variant_colors, debug)  # plot alt guides
    reflegend = _plot_reference_elements(ax, topdf, score)  # plot ref guides
    _configure_x_axis(ax, topdf, cgids_present)  # style x axis
    _configure_y_axis(ax, score)  # style y axis and title
    _add_legend(reflegend)  # add legend
    _apply_plot_styling(ax)  # style plot
    _adjust_layout(f)
    _save_plot(prefix, score, outdir)  # save and close
    plt.close()


def _draw_delta_dotplot(
    report: pd.DataFrame,
    region: Region,
    candidate_guides: List[str],
    score: str,
    guidelen: int,
    prefix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Generates and saves a delta dot plot for guide efficiency scores in a
    genomic region.

    This function computes the delta scores for guides, creates a dot plot
    visualization, and saves the resulting figure to the output directory.

    Args:
        report: DataFrame containing guide and score information.
        region: Region object representing the genomic region.
        candidate_guides: List of candidate guide sequences.
        score: The score column to visualize.
        guidelen: Length of the guide sequences.
        prefix: Prefix for output file names.
        outdir: Directory where plots will be saved.
        verbosity: Verbosity level for logging.
        debug: Boolean flag for debug mode.

    Returns:
        None
    """
    print_verbosity(
        f"Computing delta dot plot for score: {score}", verbosity, VERBOSITYLVL[3]
    )
    start = time()
    cgids = _check_candidate_guides(
        region, initialize_candidate_guides(candidate_guides, guidelen, debug), report
    )  # check if candidate guides have been provided
    score_table = _compute_delta_table(report, cgids, score)  # compute scores table
    _dotplot_delta(score_table, score, cgids, prefix, outdir, debug)  # plot delta
    print_verbosity(
        f"Delta plot computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )


def deltaplot_guides_score(
    report: pd.DataFrame,
    region: Region,
    candidate_guides: List[str],
    guidelen: int,
    prefix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Generates and saves delta dot plots for guide efficiency scores across
    multiple scoring methods.

    This function iterates over all defined score columns, generating a delta
    dot plot for each available score in the report.

    Args:
        report: DataFrame containing guide and score information.
        region: Region object representing the genomic region.
        candidate_guides: List of candidate guide sequences.
        guidelen: Length of the guide sequences.
        prefix: Prefix for output file names.
        outdir: Directory where plots will be saved.
        verbosity: Verbosity level for logging.
        debug: Boolean flag for debug mode.

    Returns:
        None
    """
    colnames = report.columns.tolist()
    for score in SCORES:  # iterate over scores to create delta plots
        if score not in colnames:  # skip delta plot
            warning(
                f"Skipping delta scores graphical report generation for score {score}",
                verbosity,
            )
            continue
        _draw_delta_dotplot(
            report,
            region,
            candidate_guides,
            score,
            guidelen,
            prefix,
            outdir,
            verbosity,
            debug,
        )


def draw_plots(
    report: pd.DataFrame,
    region: Region,
    candidate_guides: List[str],
    guidelen: int,
    prefix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    """Generates and saves graphical plots for guide types and efficiency scores.

    This function creates a pie chart of guide types and delta dot plots for guide
    efficiency scores for a given region, saving the resulting figures to the
    output directory.

    Args:
        report: DataFrame containing guide and score information.
        region: Region object representing the genomic region.
        candidate_guides: List of candidate guide sequences.
        guidelen: Length of the guide sequences.
        prefix: Prefix for output file names.
        outdir: Directory where plots will be saved.
        verbosity: Verbosity level for logging.
        debug: Boolean flag for debug mode.

    Returns:
        None
    """
    # draw guide types piechart
    piechart_guides_type(report.copy(), region, prefix, outdir, verbosity, debug)
    # draw delta plots on efficiency scores
    deltaplot_guides_score(
        report.copy(),
        region,
        candidate_guides,
        guidelen,
        prefix,
        outdir,
        verbosity,
        debug,
    )


def compute_graphical_reports(
    reports: Dict[Region, str], args: CrisprHawkSearchInputArgs
) -> None:
    """Generates graphical reports for each genomic region and saves them to disk.

    This function processes report files for multiple regions, generates graphical
    plots for each, and saves the resulting figures in the specified output directory.

    Args:
        reports: Dictionary mapping Region objects to report file paths.
        args: CrisprHawkSearchInputArgs object containing configuration parameters.

    Returns:
        None
    """
    # create figures folder in output directory
    outdir_gr = create_figures_dir(args.outdir)
    # start graphical reports computation
    print_verbosity("Computing graphical reports", args.verbosity, VERBOSITYLVL[1])
    start = time()
    for region, report in reports.items():
        # format region name as plots" name prefix
        prefix = format_region_prefix(region)
        report_df = pd.read_csv(report, sep="\t")  # load report tsv
        draw_plots(
            report_df,
            region,
            args.candidate_guides,
            args.guidelen,
            prefix,
            outdir_gr,
            args.verbosity,
            args.debug,
        )  # draw graphical reports
    print_verbosity(
        f"Graphical reports computed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
