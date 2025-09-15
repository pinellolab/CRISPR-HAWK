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

GUIDETYPES = {0: "Reference Guides", 1: "Spacer+PAM Alternative Guides", 2: "Spacer Alternative Guides", 3: "PAM Alternative Guides"}
SCORES = ["score_azimuth", "score_rs3", "score_deepcpf1", "score_cfdon", "score_elevationon"]


def _format_region(region: Region) -> str:
    return f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}"

def compute_guide_id(chrom: str, start: int, stop: int, strand: str, sgrna: str, pam: str) -> str:
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
    exception_handler(CrisprHawkGraphicalReportsError, f"Unknown guide type for grna {sgrna}", os.EX_DATAERR, debug)

def _count_guide_type(guide_types: List[int]) -> Dict[str, int]:
    types_data = {label: 0 for _, label in GUIDETYPES.items()}
    for gt in guide_types:  # count number of guides for each type
        types_data[GUIDETYPES[gt]] += 1
    return types_data     

def _draw_piechart(data: Dict[str, int], region_format: str, prefix: str, outdir: str) -> None:
    labels = list(data.keys())  # pie chart labels
    values = list(data.values())  # pie chart data
    colors = ["#5f8dd3ff", "#0055d4ff", "#ff6600ff", "#ffcc00ff"]
    explode = (0, 0, 0.05, 0.1)
    f, ax = plt.subplots(1, 1, figsize=(10, 10))
    wedges, texts, autotexts = ax.pie(values, explode=explode, colors=colors, autopct="%.2f", shadow=False, startangle=140, textprops={"fontsize": 14}, pctdistance=1.1)
    ax.set_title(f"Guide Types (region: {region_format})", fontsize=20)
    plt.legend(labels, loc=(0.8, 0), prop={"size": 16})
    plt.axis("equal")
    plt.tight_layout()
    piechart_fname = os.path.join(outdir, f"{prefix}_guides_type.png")
    plt.savefig(piechart_fname, format="png", dpi=300)


def piechart_guides_type(report: pd.DataFrame, region: Region, prefix: str, outdir: str, verbosity: int, debug: bool) -> None:
    print_verbosity("Computing guide type pie chart", verbosity, VERBOSITYLVL[3])
    start = time()
    # compute guide ids and drop non unique sites
    report["guide_id"] = report.apply(lambda x: compute_guide_id(x["chr"], x["start"], x["stop"], x["strand"], x["sgRNA_sequence"], x["pam"]), axis=1)
    report = report.drop_duplicates(subset="guide_id")
    # assign guide type 
    report["guide_type"] = report.apply(lambda x: assign_guide_type(x["origin"], x["sgRNA_sequence"], x["pam"], debug), axis=1)
    # draw pie chart for guide types
    _draw_piechart(_count_guide_type(report["guide_type"].tolist()), str(region.coordinates), prefix, outdir)
    print_verbosity(f"Guide type pie chart computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])



def compute_guide_coord(chrom: str, start: int, stop: int, strand: str) -> str:
    return f"{chrom}_{start}_{stop}_{strand}"


def _add_guide_ids(df: pd.DataFrame, debug: bool) -> pd.DataFrame:
    try:
        df["guide_id"] = df.apply(
            lambda x: compute_guide_coord(x.iloc[0], x.iloc[1], x.iloc[2], x.iloc[6]), 
            axis=1
        )
        return df
    except (IndexError, KeyError) as e:
        exception_handler(CrisprHawkGraphicalReportsError, f"An error occurred while generating guide IDs: {e}", os.EX_DATAERR, debug)


def _compute_group_delta(group: pd.DataFrame, score: str) -> pd.DataFrame:
    reference_row = group[group["origin"] == "ref"]
    if reference_row.empty:
        return group  # no reference, leave delta to 0
    reference_score = reference_row[score].values[0]
    group["delta"] = 0.0 if len(group) == 1 else group[score] - reference_score
    return group


def calculate_deltas(report: pd.DataFrame, score: str) -> pd.DataFrame:
    report["delta"] = 0.0  # initialize delta scores
    report = report.groupby("guide_id", group_keys=False).apply(lambda group: _compute_group_delta(group, score))
    return report


def compute_delta_score_dotplot(report: pd.DataFrame, score: str, prefix: str, outdir: str) -> None:
    report["guide_id"] = report.apply(lambda x: compute_guide_coord(x[0], x[1], x[2], x[6]), axis=1)
    report = calculate_deltas(report, score)




def compute_score_table(report: pd.DataFrame, score: str, debug: bool) -> pd.DataFrame:
    # validate input data
    required_columns = ["origin", "sgRNA_sequence", "pam", score, "variant_id"]
    if any(col not in report.columns for col in required_columns):
        exception_handler(CrisprHawkGraphicalReportsError, "Missing required columns from input report", os.EX_DATAERR, debug)
    if report.empty:
        exception_handler(CrisprHawkGraphicalReportsError, "Input report is empty", os.EX_DATAERR, debug)
    df = report.copy()  # make a copy to avoid modifying the original DataFrame
    df = _add_guide_ids(df, debug)   # generate guide IDs
    
    # Handle sg1617 special case
    include_sg1617, sg1617_guide_id = _check_sg1617_inclusion(df)
    
    # Calculate delta scores
    df = calculate_deltas(df, score)
    
    # Get appropriate samples column name
    samples_col = _get_samples_column_name(dataset_type)
    
    # Process guides and build guide data dictionary
    guide_rows = _build_guide_data(df, score, samples_col, include_sg1617, sg1617_guide_id)
    
    if not guide_rows:
        raise ValueError("No valid guides found after processing")
    
    # Rank guides by worst delta and select top N
    final_guides = _select_top_guides_by_delta(
        guide_rows, top_n, include_sg1617, sg1617_guide_id
    )
    
    # Build wide-format output table
    return _build_output_table(final_guides, guide_rows, samples_col)





def _check_sg1617_inclusion(df: pd.DataFrame) -> Tuple[bool, Optional[str]]:
    """Check if sg1617 guide should be included and return its ID."""
    include_sg1617 = 'chr2' in df.iloc[:, 0].unique()
    sg1617_guide_id = None
    if include_sg1617:
        sg1617_guide_id = compute_guide_coord('chr2', 60495261, 60495284, '-')
    return include_sg1617, sg1617_guide_id


def _get_samples_column_name(dataset_type: str) -> str:
    """Get appropriate samples column name based on dataset type."""
    return 'n_samples' if dataset_type == 'gnomAD' else 'samples'


def _build_guide_data(
    df: pd.DataFrame, 
    score_col: str, 
    samples_col: str,
    include_sg1617: bool,
    sg1617_guide_id: Optional[str]
) -> Dict[str, Dict[str, Any]]:
    """
    Build dictionary containing reference and alternative data for each guide.
    
    Args:
        df: DataFrame with guide data
        score_col: Name of score column
        samples_col: Name of samples column
        include_sg1617: Whether to include sg1617 guide
        sg1617_guide_id: Guide ID for sg1617
        
    Returns:
        Dictionary mapping guide IDs to their data
    """
    grouped = df.groupby('guide_id', sort=False)
    guide_rows = {}
    
    for guide_id, group in grouped:
        guide_data = _process_single_guide(
            group, score_col, samples_col, guide_id, 
            include_sg1617, sg1617_guide_id
        )
        
        if guide_data is not None:
            guide_rows[guide_id] = guide_data
    
    return guide_rows


def _process_single_guide(
    group: pd.DataFrame,
    score_col: str,
    samples_col: str,
    guide_id: str,
    include_sg1617: bool,
    sg1617_guide_id: Optional[str]
) -> Optional[Dict[str, Any]]:
    """
    Process a single guide group and return its data dictionary.
    
    Args:
        group: DataFrame group for a single guide
        score_col: Name of score column
        samples_col: Name of samples column
        guide_id: Current guide ID
        include_sg1617: Whether to include sg1617 guide
        sg1617_guide_id: Guide ID for sg1617
        
    Returns:
        Dictionary with guide data or None if invalid
    """
    ref = group[group['origin'] == 'ref']
    alts = group[group['origin'] == 'alt']
    
    if ref.empty:
        return None
    
    # Extract reference information
    ref_score = ref[score_col].values[0]
    ref_seq = ref['sgRNA_sequence'].values[0]
    ref_pam = ref['pam'].values[0]
    
    # Find valid alternatives (scores lower than reference)
    valid_alts = alts[alts[score_col] < ref_score]
    
    # Handle sg1617 special case - always include even without valid alts
    if (include_sg1617 and guide_id == sg1617_guide_id and valid_alts.empty):
        return {
            'ref_sgRNA': ref_seq,
            'pam': ref_pam,
            'ref_score': ref_score,
            'alts': []
        }
    
    # Build alternatives list
    alt_list = []
    if not valid_alts.empty:
        alt_list = _build_alternatives_list(valid_alts, score_col, samples_col)
    
    return {
        'ref_sgRNA': ref_seq,
        'pam': ref_pam,
        'ref_score': ref_score,
        'alts': alt_list
    }


def _build_alternatives_list(
    valid_alts: pd.DataFrame, 
    score_col: str, 
    samples_col: str
) -> List[Dict[str, Any]]:
    """Build list of alternative sequences data."""
    alt_list = []
    for _, alt in valid_alts.iterrows():
        alt_list.append({
            'alt_sgRNA': alt['sgRNA_sequence'],
            'pam': alt['pam'],
            'alt_score': alt[score_col],
            'delta': alt['delta'],
            samples_col: alt[samples_col],
            'variant_id': alt['variant_id'],
        })
    return alt_list


def _select_top_guides_by_delta(
    guide_rows: Dict[str, Dict[str, Any]], 
    top_n: int,
    include_sg1617: bool,
    sg1617_guide_id: Optional[str]
) -> pd.DataFrame:
    """
    Select top N guides based on worst delta scores.
    
    Args:
        guide_rows: Dictionary of guide data
        top_n: Number of guides to select
        include_sg1617: Whether to include sg1617 guide
        sg1617_guide_id: Guide ID for sg1617
        
    Returns:
        DataFrame with selected guides and their ranks
    """
    # Calculate worst delta for each guide
    worst_deltas = []
    for guide_id, data in guide_rows.items():
        if not data['alts']:
            worst_delta = 0.0
        else:
            worst_delta = min(alt['delta'] for alt in data['alts'])
        worst_deltas.append((guide_id, worst_delta))
    
    # Create DataFrame and sort by delta
    worst_df = pd.DataFrame(
        worst_deltas, 
        columns=['guide_id', 'delta']
    ).sort_values('delta').reset_index(drop=True)
    
    # Handle sg1617 special selection logic
    if include_sg1617 and sg1617_guide_id:
        final_guides = _handle_sg1617_selection(worst_df, top_n, sg1617_guide_id)
    else:
        final_guides = worst_df.head(top_n).copy()
        final_guides['Rank'] = final_guides.index + 1
    
    return final_guides


def _handle_sg1617_selection(
    worst_df: pd.DataFrame, 
    top_n: int, 
    sg1617_guide_id: str
) -> pd.DataFrame:
    """Handle special selection logic for sg1617 guide."""
    sg1617_mask = worst_df['guide_id'] == sg1617_guide_id
    
    if sg1617_mask.any():
        # Include sg1617 and top (n-1) others
        sg1617_row = worst_df[sg1617_mask].iloc[0]
        others = worst_df[~sg1617_mask].head(top_n - 1)
        final_guides = pd.concat([
            pd.DataFrame([sg1617_row]), others
        ], ignore_index=True)
        final_guides['Rank'] = final_guides.index + 1
    else:
        # sg1617 not found, just take top n
        final_guides = worst_df.head(top_n).copy()
        final_guides['Rank'] = final_guides.index + 1
    
    return final_guides


def _build_output_table(
    final_guides: pd.DataFrame, 
    guide_rows: Dict[str, Dict[str, Any]], 
    samples_col: str
) -> pd.DataFrame:
    """
    Build wide-format output table with guide data.
    
    Args:
        final_guides: DataFrame with selected guides
        guide_rows: Dictionary of guide data
        samples_col: Name of samples column
        
    Returns:
        Wide-format DataFrame ready for visualization
    """
    # Determine maximum number of alternatives across all selected guides
    max_alts = max(
        len(guide_rows[guide_id]['alts']) 
        for guide_id in final_guides['guide_id']
    ) if not final_guides.empty else 0
    
    rows = []
    for _, row in final_guides.iterrows():
        guide_id = row['guide_id']
        data = guide_rows[guide_id]
        
        # Build output row starting with basic guide info
        out_row = {
            'guide_id': guide_id,
            'Rank': row['Rank'],
            'ref_sgRNA': data['ref_sgRNA'],
            'pam': data['pam'],
            'ref_score': data['ref_score'],
        }
        
        # Add alternative columns
        _add_alternative_columns(out_row, data['alts'], max_alts, samples_col)
        rows.append(out_row)
    
    return pd.DataFrame(rows)


def _add_alternative_columns(
    out_row: Dict[str, Any], 
    alts: List[Dict[str, Any]], 
    max_alts: int, 
    samples_col: str
) -> None:
    """Add alternative sequence columns to output row."""
    for i in range(max_alts):
        alt_prefix = f'alt{i+1}'
        
        if i < len(alts):
            # Add data for existing alternative
            alt = alts[i]
            out_row.update({
                f'{alt_prefix}_sgRNA': alt['alt_sgRNA'],
                f'{alt_prefix}_pam': alt['pam'],
                f'{alt_prefix}_score': alt['alt_score'],
                f'{alt_prefix}_delta': alt['delta'],
                f'{alt_prefix}_{samples_col}': alt[samples_col],
                f'{alt_prefix}_variant_id': alt['variant_id']
            })
        else:
            # Fill with NaN for missing alternatives
            out_row.update({
                f'{alt_prefix}_sgRNA': np.nan,
                f'{alt_prefix}_pam': np.nan,
                f'{alt_prefix}_score': np.nan,
                f'{alt_prefix}_delta': np.nan,
                f'{alt_prefix}_{samples_col}': np.nan,
                f'{alt_prefix}_variant_id': np.nan
            })








def dotplot_delta(df, gene, score_col, top_n, prefix, outdir, dataset_type='samples'):
    top_df = df.sort_values('Rank').head(top_n).copy()
    plt.figure(figsize=(25, 8))

    top_df['rank_chr_start'] = top_df.apply(
        lambda row: f"Rank {row['Rank']}, {row['guide_id'].split('_')[0]}:{row['guide_id'].split('_')[1]}", axis=1
    )
    if dataset_type == 'gnomAD':
        s = 'n_samples'
    else:
        s = 'samples'

    # Collect variant keys across all alt columns
    alt_cols = [c for c in df.columns if c.startswith("alt") and c.endswith(f"_{s}")]
    variant_keys = []
    for _, row in top_df.iterrows():
        for col in alt_cols:
            if pd.notna(row[col]) and row[col] != 'REF':
                variant_keys.append(row['rank_chr_start'])
                break
    variant_keys = list(set(variant_keys))  # unique

    # Colors for variants
    base_cmaps = ['Purples', 'Blues', 'Greens', 'Oranges', 'Reds', 'PuRd', 'RdPu', 'BuPu',
                  'GnBu', 'PuBuGn', 'BuGn', 'Spectral', 'coolwarm']
    palette = [sns.color_palette(cmap, 9)[7] for cmap in base_cmaps]

    if len(variant_keys) > len(palette):
        extra_colors = []
        for i in range(6, 9):
            for cmap in base_cmaps:
                extra_colors.append(sns.color_palette(cmap, 9)[i])
        palette += extra_colors[:len(variant_keys) - len(palette)]

    random.shuffle(palette)
    variant_colors = {k: palette[i] for i, k in enumerate(variant_keys)}

    sg1617_guide_id = compute_guide_coord('chr2', 60495261, 60495284, '-')

    # Plotting
    for _, row in top_df.iterrows():
        x = row['Rank']

        # Reference point
        if row['guide_id'] == sg1617_guide_id:
            plt.scatter(
                x, row['ref_score'],
                color='gray', s=300, alpha=1,
                edgecolors='black', linewidth=1.2,
                marker='*', label='sg1617 (Baseline)' if x == 0 else ""
            )
        else:
            plt.scatter(
                x, row['ref_score'],
                color='black', s=250, alpha=0.6,
                edgecolors='white', linewidth=0.5,
                label='Reference' if x == 1 else ""
            )

        # Loop over all alts for this guide
        for i in range(1, 1000):  # upper bound
            score_col_alt = f"alt{i}_score"
            samples_col_alt = f"alt{i}_{s}"
            if score_col_alt not in row or samples_col_alt not in row:
                break
            if pd.isna(row[score_col_alt]):
                continue

            color = variant_colors.get(row['rank_chr_start'], 'gray')
            
            # Handle sample count based on dataset type
            if dataset_type == 'gnomAD':
                n_samples = int(row[samples_col_alt]) if pd.notna(row[samples_col_alt]) else 1
            else:
                n_samples = len(str(row[samples_col_alt]).split(',')) if row[samples_col_alt] != 'REF' else 1
            
            size = 150 * np.sqrt(n_samples)
            marker = 'D' if n_samples == 1 else 'o'

            plt.scatter(
                x, row[score_col_alt],
                color=color, s=size, alpha=0.6,
                edgecolors='white', linewidth=0.5,
                marker=marker
            )

    plt.xticks(top_df['Rank'], labels=top_df['ref_sgRNA'], rotation=90, fontsize=14)
    plt.xlabel('Rank', fontsize=18)
    ylabel = 'On Target Score for Mismatching Guides' if score_col == 'score_cfdon' else score_col
    plt.ylabel(ylabel, fontsize=18)
    plt.title(f'Variant-Associated Drop in Editing Efficiency - {gene}', fontsize=22)

    # Legends
    ref_patch = plt.Line2D([0], [0], marker='o', color='w', label='Reference',
                           markerfacecolor='black', markersize=np.sqrt(250), alpha=0.6)
    sg1617_patch = plt.Line2D([0], [0], marker='*', color='w', label='sg1617',
                              markerfacecolor='gray', markersize=np.sqrt(300),
                              alpha=1, markeredgecolor='black', markeredgewidth=1.2)

    legend_items = sorted(
        variant_colors.items(),
        key=lambda x: int(x[0].split(',')[0].replace('Rank ', '').strip())
    )
    representative_n_samples = 30
    variant_size = 150 * np.sqrt(representative_n_samples)
    variant_patches = [
        plt.Line2D([0], [0], marker='o', color='w',
                   label=label, markerfacecolor=color,
                   markersize=np.sqrt(variant_size) * 0.8, alpha=0.6)
        for label, color in legend_items
    ]

    # Conditional legend handles based on sg1617 presence
    legend_handles = [ref_patch]
    if sg1617_guide_id in df['guide_id'].values:
        legend_handles.append(sg1617_patch)
    legend_handles += variant_patches

    main_legend = plt.legend(
        handles=legend_handles,
        frameon=False, bbox_to_anchor=(0.35, -0.55),
        loc='upper center',
        title='Dot Color = Guide (Rank, Chrom:Start)',
        title_fontsize=15, fontsize=12, ncol=6,
        handletextpad=1.5, labelspacing=2, borderpad=1.2
    )
    plt.gca().add_artist(main_legend)

    # Size legend
    legend_labels = ['1 sample', '2-20 samples', '21-50 samples',
                     '51-100 samples', '101-200 samples', '> 200 samples']
    sample_counts = [1, 10, 35, 75, 150, 300]
    markers = ['D', 'o', 'o', 'o', 'o', 'o']
    scaled_sizes = [150 * np.sqrt(n) for n in sample_counts]

    size_legend_handles = [
        plt.Line2D([0], [0], marker=marker, color='w', label=label,
                   markerfacecolor='gray', markersize=np.sqrt(size), alpha=0.6)
        for label, size, marker in zip(legend_labels, scaled_sizes, markers)
    ]

    plt.legend(
        handles=size_legend_handles,
        frameon=False, bbox_to_anchor=(0.95, -0.55),
        loc='upper center',
        title='Dot Size = # Samples',
        title_fontsize=15, fontsize=12, ncol=2,
        handletextpad=2.5, labelspacing=3, borderpad=1.2
    )

    if score_col == 'score_deepcpf1':
        plt.ylim(-0.2, 102)
    else:
        plt.ylim(-0.2, 1.2)
    
    plt.grid(True, alpha=0.3, linestyle='--')
    sns.despine()
    dotplot_fname = os.path.join(outdir, f"{prefix}_{score_col}_delta.png")
    plt.savefig(dotplot_fname, format="png", dpi=300)


    


def compute_graphical_reports(reports: Dict[Region, str], outdir: str, verbosity: int, debug: bool) -> None:
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
                warning(f"Skipping score variation graphical report generation for score {score}", VERBOSITYLVL[2])
                continue
            print_verbosity(f"Computing delta dot plot for score: {score}", verbosity, VERBOSITYLVL[3])
            start = time()
            dotplot_delta(prepare_data_by_delta(report_df, score, 25, "1000G"), "BCL11A", score, 25, prefix, outdir_gr)
            print_verbosity(f"Delta plot computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    print_verbosity(
        f"Graphical reports computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )





