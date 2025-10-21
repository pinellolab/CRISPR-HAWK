""" """

from .crisprhawk_error import CrisprHawkCandidateGuideError
from .exception_handlers import exception_handler
from .coordinate import Coordinate
from .utils import warning, print_verbosity, CANDIDATEGUIDESREPORTPREFIX, VERBOSITYLVL
from .reports import REPORTCOLS
from .region import Region
from .pam import PAM

from typing import List, Dict
from time import time

import pandas as pd

import os

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
        exception_handler(CrisprHawkCandidateGuideError, f"Forbidden candidate guide coordinate ({coordinate})", os.EX_DATAERR, debug, e)


def initialize_candidate_guides(candidate_guides: List[str], guidelen: int, debug: bool) -> List[CandidateGuide]:
    return [CandidateGuide(g, guidelen, debug) for g in candidate_guides]

def initialize_region_reports(report_fnames: Dict[Region, str]) -> Dict[Coordinate, str]:
    return {region.coordinates: fname for region, fname in report_fnames.items()}

def _subset_region_report(report: pd.DataFrame, cg: CandidateGuide, debug: bool) -> pd.DataFrame:
    report_sub = report[report[REPORTCOLS[1]] == cg.position]  # only alternative to candidate
    if report_sub.empty:  # empty?
        exception_handler(CrisprHawkCandidateGuideError, f"Candidate guide {cg} not found. Is the candidate guide correct?", os.EX_DATAERR, debug)
    report_sub.reset_index(drop=True, inplace=True)  # reset index
    return report_sub

def _store_region_report_subset(report_sub: pd.DataFrame, cg: CandidateGuide, pam: PAM, guidelen: int, outdir: str) -> str:
    report_sub_fname = os.path.join(outdir, f"{CANDIDATEGUIDESREPORTPREFIX}__{cg.contig}_{cg.position}_{pam}_{guidelen}.tsv")
    report_sub.to_csv(report_sub_fname, sep="\t", index=False, na_rep="NA")  # store sub report
    return report_sub_fname


def subset_reports(candidate_guides: List[CandidateGuide], region_reports: Dict[Coordinate, str], pam: PAM, guidelen: int, outdir: str, debug: bool) -> Dict[CandidateGuide, str]:
    cg_reports = {}  # candidate guides sub reports 
    for region, report_fname in region_reports.items():
        report = pd.read_csv(report_fname, sep="\t")  # load region grnas report
        for cg in candidate_guides:  # retrieve candidate guides in current region
            if region.contains(cg.coordinate):
                report_sub = _subset_region_report(report, cg, debug)
                if not report_sub.empty:
                    cg_reports[cg] = _store_region_report_subset(report_sub, cg, pam, guidelen, outdir)
    return cg_reports

    
def candidate_guides_analysis(candidate_guides_coords: List[str], reports: Dict[Region, str], pam: PAM, guidelen: int, outdir: str, verbosity: int, debug: bool) -> None:
    print_verbosity("Analyzing candidate guides", verbosity, VERBOSITYLVL[1])
    start = time()
    candidate_guides = initialize_candidate_guides(candidate_guides_coords, guidelen, debug)  # initialize canididate guide objects
    region_reports = initialize_region_reports(reports)  # create report map
    # subset reports for candidate guides
    cg_reports = subset_reports(candidate_guides, region_reports, pam, guidelen, outdir, debug)
    print(cg_reports)
    print_verbosity(
        f"Candidate guides analysis completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )

