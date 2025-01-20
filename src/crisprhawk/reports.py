"""
"""

from bedfile import Region, IndelRegion
from search_guides import match
from crisprhawk_error import CrisprHawkBitsetError, CrisprHawkGuidesReportError
from exception_handlers import exception_handler
from sequences import PAM
from enrichment import INDELTYPES
from utils import GUIDESREPORTPREFIX, IUPACTABLE, IUPAC, reverse_complement

from typing import Tuple, List, Union, Dict, Any

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


def compute_pam_class(pam: PAM) -> str:
    # retrieve a string representing the input pam class
    # e.g. NGG -> [ACGT]GG
    return "".join([nt if nt in IUPAC[:4] else f"[{IUPACTABLE[nt]}]" for nt in pam.pam])


def adjust_coordinates(
    position: int,
    guidelen: int,
    pamlen: int,
    region: Union[Region, IndelRegion],
    strand: str,
) -> Tuple[int, int]:
    # adjust start position (matching position is placed on pam start)
    start = (
        (position - guidelen) + region.start + 1
        if strand == "+"
        else position + region.start + 1
    )
    stop = start + guidelen + pamlen  # adjust stop position
    if isinstance(region, IndelRegion):  # adjust positions according to indel
        stop = (
            stop - region.indel_len
            if region.indel_type == INDELTYPES[0]
            else stop + region.indel_len
        )
    return start, stop


def compute_guide_origin(samples: str) -> str:
    # compute whether the guide came from reference or alternative genomes
    origin = "ref" if samples == "NA" else "alt"
    return origin


def update_report_fields(
    report: Dict[str, List[Any]],
    start: int,
    stop: int,
    guideseq: str,
    pamguide: str,
    pam: PAM,
    strand: str,
    samples: str,
    variants: str,
    region: Union[Region, IndelRegion],
    guidelen: int,
) -> Dict[str, List[Any]]:
    # update report fields
    report[REPORTCOLS[1]].append(start)  # start position
    report[REPORTCOLS[2]].append(stop)  # stop position
    report[REPORTCOLS[3]].append(guideseq)  # guide
    report[REPORTCOLS[4]].append(pamguide)  # pam guide
    # compute extended pam class for the input pam
    report[REPORTCOLS[5]].append(compute_pam_class(pam))
    report[REPORTCOLS[6]].append(strand)  # strand orientation
    report[REPORTCOLS[7]].append(compute_guide_origin(samples))  # genome of origin
    report[REPORTCOLS[8]].append(samples)  # samples list
    report[REPORTCOLS[9]].append(variants)  # variant ids
    report[REPORTCOLS[10]].append(region.format(pad=guidelen))  # target region
    return report


def process_data(
    region: Union[Region, IndelRegion],
    guides: List[List[str]],
    matches: List[int],
    samples: List[str],
    variants: List[str],
    pam: PAM,
    strand: str,
    right: bool,
    guidelen: int,
    debug: bool,
) -> pd.DataFrame:
    report = {cname: [] for cname in REPORTCOLS}  # initialize report dictionary
    pamlen = len(pam)  # pam length used to extract pam sequence
    for i, guide in enumerate(guides):  # add guides data to report
        g = reverse_complement(guide, debug) if strand == "-" else guide
        pamguide = g[:pamlen] if right else g[-pamlen:]
        guideseq = g[pamlen:] if right else g[:-pamlen]
        if keepguide(pam, pamguide.upper(), debug):
            report[REPORTCOLS[0]].append(region.contig)  # chromosome
            # compute start and stop positions wrt region
            # TODO: handle starting on indel
            start, stop = adjust_coordinates(
                matches[i], guidelen, pamlen, region, strand
            )
            # update report fields
            report = update_report_fields(
                report,
                start,
                stop,
                guideseq,
                pamguide,
                pam,
                strand,
                samples[i],
                variants[i],
                region,
                guidelen,
            )
    return pd.DataFrame(report)


def construct_report(
    region: Union[Region, IndelRegion],
    guides: Tuple[List[List[str]], List[List[str]]],
    matches: Tuple[List[int], List[int]],
    samples: Tuple[List[str], List[str]],
    variants: Tuple[List[str], List[str]],
    pam: PAM,
    right: bool,
    guidelen: int,
    debug: bool,
) -> pd.DataFrame:
    # compute single reports on forward strand guides and reverse strand guides
    report = pd.concat(
        [
            process_data(  # forward strand
                region,
                guides[0],
                matches[0],
                samples[0],
                variants[0],
                pam,
                "+",
                right,
                guidelen,
                debug,
            ),
            process_data(  # reverse strand
                region,
                guides[1],
                matches[1],
                samples[1],
                variants[1],
                pam,
                "-",
                right,
                guidelen,
                debug,
            ),
        ]
    )
    return report


def split_region_reports(
    guides_raw: Dict[
        Union[Region, IndelRegion],
        Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]],
    ],
    pam: PAM,
    right: bool,
    guidelen: int,
    debug: bool,
) -> Tuple[Dict[Region, pd.DataFrame], Dict[Region, pd.DataFrame]]:
    # dictionaries storing reports for query and indel regions
    reports_regions, reports_indels = {}, {}
    for region, (positions, guides, samples, variants) in guides_raw.items():
        # construct individual region report
        report = construct_report(
            region, guides, positions, samples, variants, pam, right, guidelen, debug
        )
        if isinstance(region, IndelRegion):  # indel region
            reports_indels[region] = report
        else:  # query region
            reports_regions[region] = report
    return reports_regions, reports_indels


def merge_reports(
    reports_regions: Dict[Region, pd.DataFrame],
    reports_indels: Dict[IndelRegion, pd.DataFrame],
) -> Dict[Region, pd.DataFrame]:
    # for each indel report search the original query region
    for region_indel, report_indel in reports_indels.items():
        for region, report in reports_regions.items():
            if region.contains(region_indel):
                reports_regions[region] = pd.concat([report, report_indel])
                break  # once found the query region, go to next region
    return reports_regions


def format_report(report: pd.DataFrame, region_name: str) -> pd.DataFrame:
    report[REPORTCOLS[10]] = region_name  # force query region report
    # reset dataframe index and sort by genomic coordinates
    report = report.reset_index(drop=True)
    report = report.sort_values([REPORTCOLS[0], REPORTCOLS[1]], ascending=True)
    # force start and stop to int values - they may be treated as float if
    # concatenated with empty dataframe (e.g. no guide found on + or - strand)
    report[REPORTCOLS[1]] = report[REPORTCOLS[1]].astype(int)
    report[REPORTCOLS[2]] = report[REPORTCOLS[2]].astype(int)
    return report


def store_reports(
    reports: Dict[Region, pd.DataFrame],
    guidelen: int,
    pam: PAM,
    outdir: str,
    debug: bool,
) -> None:
    try:
        for region, report in reports.items():
            region_name = region.format(pad=guidelen)  # query region name
            # define report fname -> example: crisprhawk_guides__chr2_1_10_NGG_20.tsv
            guidesreport = os.path.join(
                outdir, f"{GUIDESREPORTPREFIX}__{region_name}_{pam}_{guidelen}.tsv"
            )
            report = format_report(report, region_name)  # format report
            report.to_csv(guidesreport, sep="\t", index=False)  # store report
    except FileNotFoundError as e:
        exception_handler(
            CrisprHawkGuidesReportError,
            f"Unable to write to {outdir}",
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


def report_guides(
    guides_raw: Dict[
        Union[Region, IndelRegion],
        Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]],
    ],
    outdir: str,
    pam: PAM,
    right: bool,
    guidelen: int,
    debug: bool,
):
    # split reports by region type
    reports_regions, reports_indels = split_region_reports(
        guides_raw, pam, right, guidelen, debug
    )
    # merge region reports referring to the same query regions
    reports = merge_reports(reports_regions, reports_indels)
    # format and store reports
    store_reports(reports, guidelen, pam, outdir, debug)
