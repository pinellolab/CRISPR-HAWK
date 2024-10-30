"""
"""

from bedfile import Region
from search_guides import match
from crisprhawk_error import CrisprHawkBitsetError, CrisprHawkGuidesReportError
from exception_handlers import exception_handler
from sequences import PAM, explode_iupac_sequence
from utils import GUIDESREPORTPREFIX, reverse_complement

from typing import Tuple, List

import pandas as pd

import os

REPORTCOLS = [
    "Chr",
    "PAM_hit",
    "PAM_seq",
    "sgRNA_sequence",
    "Strand",
    "sgRNA_start",
    "sgRNA_stop",
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


def construct_report(
    region: Region,
    guides: List[List[str]],
    matches: List[int],
    pam: PAM,
    strand: str,
    right: bool,
    guidelen: int,
    debug: bool,
) -> pd.DataFrame:
    report = {cname: [] for cname in REPORTCOLS}  # initialize report dictionary
    pamlen = len(pam)  # pam length used to extract pam sequence
    for i, guide in enumerate(guides):  # add guides data to report
        for g in explode_iupac_sequence(guide, debug):
            g = reverse_complement(g, debug) if strand == "-" else g
            pamguide = g[:pamlen] if right else g[-pamlen:]
            guide = g[pamlen:] if right else g[:-pamlen]
            if keepguide(pam, pamguide, debug):
                report[REPORTCOLS[0]].append(region.contig)  # chromosome
                report[REPORTCOLS[1]].append(pamguide)  # pam hit
                report[REPORTCOLS[2]].append(pam.pam)  # query pam
                report[REPORTCOLS[3]].append(g)  # query pam
                report[REPORTCOLS[4]].append(strand)  # strand orientation
                # compute start and stop positions wrt region
                start = (matches[i] - guidelen) + region.start
                stop = start + guidelen + pamlen
                report[REPORTCOLS[5]].append(start)  # start position
                report[REPORTCOLS[6]].append(stop)  # stop position
    return pd.DataFrame(report)


def report_guides(
    outdir: str,
    region: Region,
    guides: Tuple[List[List[str]], List[List[str]]],
    matches: Tuple[List[int], List[int]],
    pam: PAM,
    right: bool,
    guidelen: int,
    debug: bool,
) -> None:
    # define report fname -> example: crisprhawk_guides__chr2_1_10_NGG_20.tsv
    guidesreport = os.path.join(
        outdir, f"{GUIDESREPORTPREFIX}__{region.format()}_{pam}_{guidelen}.tsv"
    )
    try:
        # compute single reports on forward strand guides and reverse strand guides
        report = pd.concat(
            [
                construct_report(
                    region, guides[0], matches[0], pam, "+", right, guidelen, debug
                ),
                construct_report(
                    region, guides[1], matches[1], pam, "-", right, guidelen, debug
                ),
            ]
        )
        # sort report by contig and position
        report = report.sort_values([REPORTCOLS[0], REPORTCOLS[5]], ascending=True)
        report.to_csv(guidesreport, sep="\t", index=False)  # save report
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
