""" """

from .crisprhawk_error import CrisprHawkGraphicalReportsError
from .region_constructor import PADDING
from .exception_handlers import exception_handler
from .utils import VERBOSITYLVL, is_lowercase, create_folder, print_verbosity, warning
from .region import Region


from typing import Dict, List, Tuple, Any, Set, Optional
from time import time

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

import random
import os


# configuration constants
GUIDETYPES = {
    0: "Reference Guides",
    1: "Spacer+PAM Alternative Guides",
    2: "Spacer Alternative Guides",
    3: "PAM Alternative Guides",
}
SCORES = [
    "score_azimuth",
    "score_rs3",
    "score_deepcpf1",
    "score_cfdon",
    "score_elevationon",
]
FIGURE_SIZE = (25, 15)
DPI = 300
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
SAMPLE_SIZE_MULTIPLIER = 150
REPRESENTATIVE_SAMPLES = 30


def _format_region(region: Region) -> str:
    return f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}"


def compute_guide_id(
    chrom: str, start: int, stop: int, strand: str, sgrna: str, pam: str
) -> str:
    return f"{chrom}_{start}_{stop}_{strand}_{sgrna}_{pam}"


def assign_guide_type(origin: str, sgrna: str, pam: str, debug: bool):
    # types: 0 -> ref; 1 -> spacer+pam alt; 2 -> spacer alt; 3 -> pam alt
    if origin == "ref":  # reference type guide
        return 0
    assert origin != "ref"  # asses alternative guide type
    sgrna_alt = is_lowercase(sgrna)  # check sgrna
    pam_alt = is_lowercase(pam)  # check pam
    if sgrna_alt and pam_alt:  # spacer+pam alt
        return 1
    if sgrna_alt and not pam_alt:  # spacer alt
        return 2
    if not sgrna_alt and pam_alt:  # pam alt
        return 3
    # this chunk of code should never be reached
    exception_handler(
        CrisprHawkGraphicalReportsError,
        f"Unknown guide type for grna {sgrna}",
        os.EX_DATAERR,
        debug,
    )


def _count_guide_type(guide_types: List[int]) -> Dict[str, int]:
    types_data = {label: 0 for _, label in GUIDETYPES.items()}
    for gt in guide_types:  # count number of guides for each type
        types_data[GUIDETYPES[gt]] += 1
    return types_data


def _draw_piechart(
    data: Dict[str, int], region_format: str, prefix: str, outdir: str
) -> None:
    labels = list(data.keys())  # pie chart labels
    values = list(data.values())  # pie chart data
    colors = ["#5f8dd3ff", "#0055d4ff", "#ff6600ff", "#ffcc00ff"]
    explode = (0, 0, 0.05, 0.1)
    f, ax = plt.subplots(1, 1, figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(
        values,
        explode=explode,
        colors=colors,
        autopct="%.2f",
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
    print_verbosity("Computing guide type pie chart", verbosity, VERBOSITYLVL[3])
    start = time()
    # compute guide ids and drop non unique sites
    report["guide_id"] = report.apply(
        lambda x: compute_guide_id(
            x["chr"], x["start"], x["stop"], x["strand"], x["sgRNA_sequence"], x["pam"]
        ),
        axis=1,
    )
    report = report.drop_duplicates(subset="guide_id")
    # assign guide type
    report["guide_type"] = report.apply(
        lambda x: assign_guide_type(x["origin"], x["sgRNA_sequence"], x["pam"], debug),
        axis=1,
    )
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


def compute_guide_coord(chrom: str, start: int, stop: int, strand: str) -> str:
    return f"{chrom}_{start}_{stop}_{strand}"


def _add_guide_ids(df: pd.DataFrame, debug: bool) -> pd.DataFrame:
    try:
        df["guide_id"] = df.apply(
            lambda x: compute_guide_coord(x.iloc[0], x.iloc[1], x.iloc[2], x.iloc[6]),
            axis=1,
        )
        return df
    except (IndexError, KeyError) as e:
        exception_handler(
            CrisprHawkGraphicalReportsError,
            f"An error occurred while generating guide IDs: {e}",
            os.EX_DATAERR,
            debug,
        )


def _compute_group_delta(group: pd.DataFrame, score: str) -> pd.DataFrame:
    reference_row = group[group["origin"] == "ref"]
    if reference_row.empty:
        return group  # no reference, leave delta to 0
    reference_score = reference_row[score].values[0]
    group["delta"] = 0.0 if len(group) == 1 else group[score] - reference_score
    return group


def calculate_deltas(report: pd.DataFrame, score: str) -> pd.DataFrame:
    report["delta"] = 0.0  # initialize delta scores
    report = report.groupby("guide_id", group_keys=False).apply(
        lambda group: _compute_group_delta(group, score)
    )
    return report


def _build_alternatives_list(
    valid_alts: pd.DataFrame, score_col: str
) -> List[Dict[str, Any]]:
    alt_list = []
    for _, alt in valid_alts.iterrows():
        alt_list.append(
            {
                "alt_sgRNA": alt["sgRNA_sequence"],
                "pam": alt["pam"],
                "alt_score": alt[score_col],
                "delta": alt["delta"],
                "samples": alt["samples"],
                "variant_id": alt["variant_id"],
            }
        )
    return alt_list


def _process_single_guide(
    group: pd.DataFrame, score_col: str
) -> Optional[Dict[str, Any]]:
    ref = group[group["origin"] == "ref"]
    alts = group[group["origin"] == "alt"]
    if ref.empty:
        return None
    # extract reference information
    ref_score = ref[score_col].values[0]
    ref_seq = ref["sgRNA_sequence"].values[0]
    ref_pam = ref["pam"].values[0]
    # find valid alternatives (scores lower than reference)
    valid_alts = alts[alts[score_col] < ref_score]
    # handle sg1617 special case - always include even without valid alts
    if valid_alts.empty:
        return {
            "ref_sgRNA": ref_seq,
            "pam": ref_pam,
            "ref_score": ref_score,
            "alts": [],
        }
    alt_list = []  # build alternatives list
    if not valid_alts.empty:
        alt_list = _build_alternatives_list(valid_alts, score_col)
    return {
        "ref_sgRNA": ref_seq,
        "pam": ref_pam,
        "ref_score": ref_score,
        "alts": alt_list,
    }


def _build_guide_data(df: pd.DataFrame, score_col: str) -> Dict[str, Dict[str, Any]]:
    grouped = df.groupby("guide_id", sort=False)
    guide_rows = {}
    for guide_id, group in grouped:
        guide_data = _process_single_guide(group, score_col)
        if guide_data is not None:
            guide_rows[guide_id] = guide_data
    return guide_rows


def _select_top_guides_by_delta(guide_rows: Dict[str, Dict[str, Any]]) -> pd.DataFrame:
    worst_deltas = []  # calculate worst delta for each guide
    for guide_id, data in guide_rows.items():
        if not data["alts"]:
            worst_delta = 0.0
        else:
            worst_delta = min(alt["delta"] for alt in data["alts"])
        worst_deltas.append((guide_id, worst_delta))
    worst_df = (
        pd.DataFrame(worst_deltas, columns=["guide_id", "delta"])
        .sort_values("delta")
        .reset_index(drop=True)
    )  # create DataFrame and sort by delta
    final_guides = worst_df.head(25).copy()
    final_guides["Rank"] = final_guides.index + 1
    return final_guides


def _add_alternative_columns(
    out_row: Dict[str, Any], alts: List[Dict[str, Any]], max_alts: int
) -> None:
    for i in range(max_alts):
        alt_prefix = f"alt{i+1}"
        if i < len(alts):  # add data for existing alternative
            alt = alts[i]
            out_row.update(
                {
                    f"{alt_prefix}_sgRNA": alt["alt_sgRNA"],
                    f"{alt_prefix}_pam": alt["pam"],
                    f"{alt_prefix}_score": alt["alt_score"],
                    f"{alt_prefix}_delta": alt["delta"],
                    f"{alt_prefix}_samples": alt["samples"],
                    f"{alt_prefix}_variant_id": alt["variant_id"],
                }
            )
        else:  # fill with NaN for missing alternatives
            out_row.update(
                {
                    f"{alt_prefix}_sgRNA": np.nan,
                    f"{alt_prefix}_pam": np.nan,
                    f"{alt_prefix}_score": np.nan,
                    f"{alt_prefix}_delta": np.nan,
                    f"{alt_prefix}_samples": np.nan,
                    f"{alt_prefix}_variant_id": np.nan,
                }
            )


def _build_output_table(
    final_guides: pd.DataFrame, guide_rows: Dict[str, Dict[str, Any]]
) -> pd.DataFrame:
    # determine maximum number of alternatives across all selected guides
    max_alts = (
        max(len(guide_rows[guide_id]["alts"]) for guide_id in final_guides["guide_id"])
        if not final_guides.empty
        else 0
    )
    rows = []
    for _, row in final_guides.iterrows():
        guide_id = row["guide_id"]
        data = guide_rows[guide_id]
        # build output row starting with basic guide info
        out_row = {
            "guide_id": guide_id,
            "Rank": row["Rank"],
            "ref_sgRNA": data["ref_sgRNA"],
            "pam": data["pam"],
            "ref_score": data["ref_score"],
        }

        # add alternative columns
        _add_alternative_columns(out_row, data["alts"], max_alts)
        rows.append(out_row)
    return pd.DataFrame(rows)


def compute_score_table(report: pd.DataFrame, score: str, debug: bool) -> pd.DataFrame:
    # validate input data
    required_columns = ["origin", "sgRNA_sequence", "pam", score, "variant_id"]
    if any(col not in report.columns for col in required_columns):
        exception_handler(
            CrisprHawkGraphicalReportsError,
            "Missing required columns from input report",
            os.EX_DATAERR,
            debug,
        )
    if report.empty:
        exception_handler(
            CrisprHawkGraphicalReportsError,
            "Input report is empty",
            os.EX_DATAERR,
            debug,
        )
    df = report.copy()  # make a copy to avoid modifying the original DataFrame
    df = _add_guide_ids(df, debug)  # generate guide IDs
    df = calculate_deltas(df, score)  # calculate delta scores
    # process guides and build guide data dictionary
    guide_rows = _build_guide_data(df, score)
    if not guide_rows:
        exception_handler(
            CrisprHawkGraphicalReportsError,
            "No valid guides found after processing",
            os.EX_DATAERR,
            debug,
        )
    # rank guides by worst delta and select top N
    final_guides = _select_top_guides_by_delta(guide_rows)

    # Build wide-format output table
    return _build_output_table(final_guides, guide_rows)


def _get_alternative_columns(df: pd.DataFrame) -> List[str]:
    return [c for c in df.columns if c.startswith("alt") and c.endswith(f"_samples")]


def _prepare_plot_data(df: pd.DataFrame) -> Dict:
    top_df = df.sort_values("Rank").head(25).copy()
    # create rank-chromosome-start labels
    top_df["rank_chr_start"] = top_df.apply(
        lambda row: (
            f"Rank {row['Rank']}, "
            f"{row['guide_id'].split('_')[0]}:{row['guide_id'].split('_')[1]}"
        ),
        axis=1,
    )
    alt_cols = _get_alternative_columns(df)  # retrieve alternative guides columns
    return {"top_df": top_df, "samples_col": "samples", "alt_cols": alt_cols}


def _create_color_palette(n_colors: int) -> List[Tuple[float, float, float]]:
    palette = [sns.color_palette(cmap, 9)[7] for cmap in BASE_CMAPS]
    if n_colors > len(palette):  # extend palette if needed
        extra_colors = []
        for i in range(6, 9):
            for cmap in BASE_CMAPS:
                extra_colors.append(sns.color_palette(cmap, 9)[i])
        palette += extra_colors[: n_colors - len(palette)]
    random.shuffle(palette)
    return palette


def _generate_variant_colors(
    top_df: pd.DataFrame,
) -> Dict[str, Tuple[float, float, float]]:
    alt_cols = _get_alternative_columns(top_df)  # retrieve alternative guides columns
    variant_keys = set()  # collect unique variant keys
    for _, row in top_df.iterrows():
        for col in alt_cols:
            if pd.notna(row[col]) and row[col] != "REF":
                variant_keys.add(row["rank_chr_start"])
                break
    variant_keys = list(variant_keys)
    palette = _create_color_palette(len(variant_keys))  # generate color palette
    return {key: palette[i] for i, key in enumerate(variant_keys)}


def _plot_reference_point(row: pd.Series, x: int) -> None:
    plt.scatter(
        x,
        row["ref_score"],
        color="black",
        s=250,
        alpha=0.6,
        edgecolors="white",
        linewidth=0.5,
        label="Reference" if x == 1 else "",
    )


def _extract_sample_count(samples_value) -> int:
    if samples_value == "REF" or pd.isna(samples_value):
        return 1
    return len(str(samples_value).split(","))


def _plot_alternative_points(
    row: pd.Series, x: int, variant_colors: Dict[str, Tuple[float, float, float]]
) -> None:
    alt_num = 1
    max_iterations = 1000  # safety limit
    while alt_num <= max_iterations:
        score_col_alt = f"alt{alt_num}_score"
        samples_col_alt = f"alt{alt_num}_samples"
        # check if columns exist
        if score_col_alt not in row.index or samples_col_alt not in row.index:
            break
        if pd.isna(row[score_col_alt]):  # skip missing scores
            alt_num += 1
            continue
        color = variant_colors.get(row["rank_chr_start"], "gray")  # get plot parameters
        n_samples = _extract_sample_count(row[samples_col_alt])
        size = SAMPLE_SIZE_MULTIPLIER * np.sqrt(n_samples)
        marker = "D" if n_samples == 1 else "o"
        plt.scatter(x, row[score_col_alt], color=color, s=size, alpha=0.6, edgecolors="white", linewidth=0.5, marker=marker)  # type: ignore
        alt_num += 1


def _plot_data_points(
    plot_data: Dict, variant_colors: Dict[str, Tuple[float, float, float]]
) -> None:
    top_df = plot_data["top_df"]
    for _, row in top_df.iterrows():
        x = row["Rank"]
        _plot_reference_point(row, x)  # plot reference point
        _plot_alternative_points(row, x, variant_colors)  # plot alternative points


def _configure_plot_appearance(
    top_df: pd.DataFrame, region_format: str, score_col: str
) -> None:
    plt.xticks(top_df["Rank"], labels=top_df["ref_sgRNA"], rotation=90, fontsize=14)
    plt.xlabel("Rank", fontsize=16)
    ylabel = (
        "On Target Score for Mismatching Guides"
        if score_col in ["score_cfdon", "score_elevationon"]
        else score_col
    )
    plt.ylabel(ylabel, fontsize=16)
    plt.title(
        f"Variant-Associated Drop in Editing Efficiency - Region: {region_format}",
        fontsize=20,
    )


def _create_variant_patches(variant_colors: Dict[str, Tuple[float, float, float]]) -> List[plt.Line2D]:  # type: ignore
    legend_items = sorted(
        variant_colors.items(),
        key=lambda x: int(x[0].split(",")[0].replace("Rank ", "").strip()),
    )
    variant_size = SAMPLE_SIZE_MULTIPLIER * np.sqrt(REPRESENTATIVE_SAMPLES)
    return [
        plt.Line2D(  # type: ignore
            [0],
            [0],
            marker="o",
            color="w",
            label=label,
            markerfacecolor=color,
            markersize=np.sqrt(variant_size) * 0.8,
            alpha=0.6,
        )
        for label, color in legend_items
    ]


def _add_main_legend(
    variant_colors: Dict[str, Tuple[float, float, float]], df: pd.DataFrame
) -> None:
    # create basic legend handles
    ref_patch = plt.Line2D([0], [0], marker="o", color="w", label="Reference", markerfacecolor="black", markersize=np.sqrt(250), alpha=0.6)  # type: ignore
    legend_handles = [ref_patch]
    variant_patches = _create_variant_patches(variant_colors)  # add variant patches
    legend_handles.extend(variant_patches)
    # create the legend
    main_legend = plt.legend(
        handles=legend_handles,
        frameon=False,
        bbox_to_anchor=(0.35, -0.55),
        loc="upper center",
        title="Dot Color = Guide (Rank, Chrom:Start)",
        title_fontsize=15,
        fontsize=12,
        ncol=6,
        handletextpad=1.5,
        labelspacing=2,
        borderpad=1.2,
    )
    plt.gca().add_artist(main_legend)


def _add_size_legend() -> None:
    legend_config = [
        ("1 sample", 1, "D"),
        ("2-20 samples", 10, "o"),
        ("21-50 samples", 35, "o"),
        ("51-100 samples", 75, "o"),
        ("101-200 samples", 150, "o"),
        ("> 200 samples", 300, "o"),
    ]
    size_legend_handles = [
        plt.Line2D(  # type: ignore
            [0],
            [0],
            marker=marker,
            color="w",
            label=label,
            markerfacecolor="gray",
            markersize=np.sqrt(SAMPLE_SIZE_MULTIPLIER * np.sqrt(sample_count)),
            alpha=0.6,
        )
        for label, sample_count, marker in legend_config
    ]
    plt.legend(
        handles=size_legend_handles,
        frameon=False,
        bbox_to_anchor=(0.95, -0.55),
        loc="upper center",
        title="Dot Size = # Samples",
        title_fontsize=15,
        fontsize=12,
        ncol=2,
        handletextpad=2.5,
        labelspacing=3,
        borderpad=1.2,
    )


def _add_legends(
    variant_colors: Dict[str, Tuple[float, float, float]], df: pd.DataFrame
) -> None:
    _add_main_legend(variant_colors, df)
    _add_size_legend()


def _apply_plot_styling(score_col: str) -> None:
    if score_col == "score_deepcpf1":  # set y-axis limits based on score type
        plt.ylim(-0.2, 102)
    else:
        plt.ylim(-0.2, 1.2)
    plt.grid(True, alpha=0.3, linestyle="--")
    sns.despine()
    plt.tight_layout()


def _save_plot(output_prefix: str, score_col: str, output_dir: str) -> str:
    filename = f"{output_prefix}_{score_col}_delta.png"
    output_path = os.path.join(output_dir, filename)
    plt.savefig(output_path, format="png", dpi=DPI, bbox_inches="tight")
    return output_path


def _draw_delta_dotplot(
    df: pd.DataFrame, region: Region, score_col: str, prefix: str, outdir: str
) -> str:
    plot_data = _prepare_plot_data(df)  # prepare data
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)  # set up the plot
    # generate colors for variants
    variant_colors = _generate_variant_colors(plot_data["top_df"])
    _plot_data_points(plot_data, variant_colors)  # plot the data points
    # Configure plot appearance
    _configure_plot_appearance(plot_data["top_df"], str(region.coordinates), score_col)
    _add_legends(variant_colors, df)  # add legends
    _apply_plot_styling(score_col)  # set axis limits and styling
    output_path = _save_plot(prefix, score_col, outdir)  # save the plot
    plt.close(fig)
    return output_path


def compute_delta_dotplot(
    report: pd.DataFrame,
    region: Region,
    score: str,
    prefix: str,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> None:
    print_verbosity(
        f"Computing delta dot plot for score: {score}", verbosity, VERBOSITYLVL[3]
    )
    start = time()
    score_table = compute_score_table(report, score, debug)  # compute scores table
    _draw_delta_dotplot(score_table, region, score, prefix, outdir)
    print_verbosity(
        f"Delta plot computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )


def compute_graphical_reports(
    reports: Dict[Region, str], outdir: str, verbosity: int, debug: bool
) -> None:
    # create figures folder in output directory
    outdir_gr = os.path.join(outdir, "figures")
    if not os.path.isdir(outdir_gr):
        outdir_gr = create_folder(os.path.join(outdir, "figures"))
    # start graphical reports computation
    print_verbosity("Computing graphical reports", verbosity, VERBOSITYLVL[1])
    start = time()
    for region, report in reports.items():
        prefix = _format_region(region)  # format region name as plots" name prefix
        report_df = pd.read_csv(report, sep="\t")
        # draw guide types piechart
        piechart_guides_type(report_df, region, prefix, outdir_gr, verbosity, debug)
        # draw scatter plot for guide score variation
        for score in SCORES:
            if score not in report_df.columns.tolist():
                warning(
                    f"Skipping delta scores graphical report generation for score {score}",
                    VERBOSITYLVL[2],
                )
                continue
            compute_delta_dotplot(
                report_df, region, score, prefix, outdir_gr, verbosity, debug
            )
    print_verbosity(
        f"Graphical reports computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
