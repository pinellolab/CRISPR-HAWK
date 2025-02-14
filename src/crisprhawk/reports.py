"""
"""

from bedfile import Region, IndelRegion
from search_guides import match
from crisprhawk_error import CrisprHawkBitsetError, CrisprHawkGuidesReportError
from exception_handlers import exception_handler
from sequences import PAM
from guide import Guide
from enrichment import INDELTYPES
from utils import GUIDESREPORTPREFIX, IUPACTABLE, IUPAC, STRAND, reverse_complement

from typing import Tuple, List, Union, Dict, Any, Set

import pandas as pd

import os


REPORTCOLS = [
    "chr",
    "start",
    "stop",
    "sgRNA_sequence",
    "pam",
    "pam_class",
    "strand",
    "origin",
    "samples",
    "variant_id",
    "target",
]


def keepguide(pam_query: PAM, pamseq: str, debug: bool) -> bool:
    assert hasattr(pam_query.pam, "_sequence_bits")  # should already be encoded
    try:
        pam = PAM(pamseq, debug)  # construct pam object
        assert len(pam_query) == len(pam)  # their length must match
        pam.encode()  # encode pam sequence
        # if true report the sequence -> useful when considering genetic variants
        # e.g. NRG -> one pam would be NAG, the other NGG -> report only the second
        return match(pam_query.bits[0], pam.bits[0])
    except ValueError as e:
        exception_handler(
            CrisprHawkBitsetError,
            "Failed PAM encoding during report generation",
            os.EX_DATAERR,
            debug,
            e,
        )


def adjust_coordinates(guide: Guide, region: Region) -> Tuple[int, int]:
    # adjust start position (matching position is placed on pam start)
    start = guide.position + region.start + 1
    stop = start + guide.guidelen + guide.pamlen  # adjust stop position
    if isinstance(region, IndelRegion):
        stop = (
            stop - region.indel_length
            if region.indel_type == INDELTYPES[0]
            else stop + region.indel_length
        )
    return start, stop


def compute_pam_class(pam: PAM) -> str:
    # retrieve a string representing the input pam class
    # e.g. NGG -> [ACGT]GG
    return "".join([nt if nt in IUPAC[:4] else f"[{IUPACTABLE[nt]}]" for nt in pam.pam])


def compute_guide_origin(samples: str) -> str:
    # compute whether the guide came from reference or alternative genomes
    return "ref" if samples == "NA" else "alt"


def compute_strand_orientation(strand: int) -> str:
    # retrieve strand orientation
    return "+" if strand == STRAND[0] else "-"  # retrieve strand orientation


def update_report_fields(
    report: Dict[str, List[Any]], region: Region, guide: Guide, pamclass: str
) -> Dict[str, List[Any]]:
    # update report fields
    # TODO: handle starting on indel
    start, stop = adjust_coordinates(guide, region)  # compute start and stop
    report[REPORTCOLS[1]].append(start)  # start position
    report[REPORTCOLS[2]].append(stop)  # stop position
    report[REPORTCOLS[3]].append(guide.guide)  # guide sequence
    report[REPORTCOLS[4]].append(guide.pam)  # pam guide
    report[REPORTCOLS[5]].append(pamclass)  # extended pam class
    # strand orientation
    report[REPORTCOLS[6]].append(compute_strand_orientation(guide.strand))
    report[REPORTCOLS[7]].append(compute_guide_origin(guide.samples))  # genome
    report[REPORTCOLS[8]].append(guide.samples)  # samples list
    report[REPORTCOLS[9]].append(guide.variants)  # variant ids
    report[REPORTCOLS[10]].append(region.format(pad=guide.guidelen, string=True))
    return report


def process_data(
    region: Region, guides: List[Guide], pam: PAM, debug: bool
) -> pd.DataFrame:
    report = {cname: [] for cname in REPORTCOLS}  # initialize report dictionary
    pamclass = compute_pam_class(pam)  # compute extended pam class
    for guide in guides:  # iterate over guides and add to report
        if guide.strand == STRAND[1]:  # negative strand -> report 3'-5' sequence
            guide.reverse_complement()  # compute guide reverse complement sequence
        if keepguide(pam, guide.pam.upper(), debug):  # valid guide pam
            report[REPORTCOLS[0]].append(region.contig)  # region contig (chrom)
            # update report with current guide data
            report = update_report_fields(report, region, guide, pamclass)
    return pd.DataFrame(report)  # build dataframe from report data


def construct_report(
    guides: Dict[Region, List[Guide]], pam: PAM, debug: bool
) -> Dict[Region, pd.DataFrame]:
    return {
        region: process_data(region, guides_list, pam, debug)
        for region, guides_list in guides.items()
    }


def format_report(report: pd.DataFrame) -> pd.DataFrame:
    # reset dataframe index and sort by genomic coordinates
    report = report.reset_index(drop=True)
    report = report.sort_values(
        [REPORTCOLS[0], REPORTCOLS[1], REPORTCOLS[2]], ascending=True
    )
    # force start and stop to int values - they may be treated as float if
    # concatenated with empty dataframe (e.g. no guide found on + or - strand)
    report[REPORTCOLS[1]] = report[REPORTCOLS[1]].astype(int)
    report[REPORTCOLS[2]] = report[REPORTCOLS[2]].astype(int)
    return report


def store_report(report: pd.DataFrame, guidesreport: str, debug: bool) -> None:
    try:
        report = format_report(report)  # format report
        report.to_csv(guidesreport, sep="\t", index=False)  # store report
    except FileNotFoundError as e:
        exception_handler(
            CrisprHawkGuidesReportError,
            f"Unable to write to {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )
    except PermissionError as e:
        exception_handler(
            CrisprHawkGuidesReportError,
            f"Permission denied to write {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )
    except Exception as e:
        exception_handler(
            CrisprHawkGuidesReportError,
            f"An unexpected error occurred while writing {guidesreport}",
            os.EX_OSERR,
            debug,
            e,
        )


def split_reports(
    reports: Dict[Union[Region, IndelRegion], pd.DataFrame]
) -> Tuple[Dict[Region, pd.DataFrame], Dict[IndelRegion, pd.DataFrame]]:
    # split reports in reports for guides containing snps and indels
    reports_snps, reports_indels = {}, {}
    for region, report in reports.items():
        if isinstance(region, IndelRegion):  # region defined for indels + snps
            reports_indels[region] = report
        else:  # region containing just snps
            reports_snps[region] = report
    return reports_snps, reports_indels


def merge_reports(
    reports: Dict[Union[Region, IndelRegion], pd.DataFrame]
) -> Dict[Region, pd.DataFrame]:
    # retrieve indel region reports
    reports_snps, reports_indels = split_reports(reports)
    reports_merged = {}  # merged reports dictionary
    for region, report_snp in reports_snps.items():
        # retrieve indel regions fully overlapped by query region
        iregions_overlapping = [
            report_indel
            for iregion, report_indel in reports_indels.items()
            if region.contains(iregion)
        ]
        reports_merged[region] = pd.concat([report_snp] + iregions_overlapping)
        reports_merged[region][REPORTCOLS[10]] = report_snp.loc[0, REPORTCOLS[10]]
    return reports_merged


def collapse_samples(samples: pd.Series) -> str:
    return (
        ",".join(sorted(set(",".join(samples).split(",")))) if not samples.empty else ""
    )


def parse_variant_ids(variant_ids: pd.Series) -> Set[str]:
    return set(variant_ids.split(",")) if variant_ids else set()


def check_variant_ids(variant_ids_list: List[str]) -> str:
    variant_sets = [parse_variant_ids(vid) for vid in variant_ids_list]
    unique_variant_ids_sets = set(tuple(sorted(vs)) for vs in variant_sets)
    assert len(unique_variant_ids_sets) == 1
    return ",".join(sorted(unique_variant_ids_sets.pop()))


def collapse_report_entries(report: pd.DataFrame) -> pd.DataFrame:
    # Define the columns to group by
    group_cols = REPORTCOLS[:5] + REPORTCOLS[6:8]
    # Group by the key columns and apply aggregation functions
    report_collapsed = report.groupby(group_cols, as_index=False).agg(
        {
            "pam_class": "first",  # Assuming pam_class is the same across entries
            "origin": "first",  # Assuming origin does not change
            "samples": collapse_samples,  # Merge sample lists
            "variant_id": check_variant_ids,  # Ensure identical variant_id
            "target": "first",  # Keep the first target entry
        }
    )
    return report_collapsed


def report_guides(
    guides: Dict[Union[Region, IndelRegion], List[Guide]],
    guidelen: int,
    pam: PAM,
    outdir: str,
    right: bool,
    debug: bool,
) -> None:
    reports = construct_report(guides, pam, debug)  # construct reports
    reports = merge_reports(reports)  # merge indel and snp reports
    for region, report in reports.items():  # store reports in output folder
        guidesreport = os.path.join(
            outdir,
            f"{GUIDESREPORTPREFIX}__{region.format(pad=guidelen)}_{pam}_{guidelen}.tsv",
        )
        report = collapse_report_entries(report)
        store_report(report, guidesreport, debug)  # write report
