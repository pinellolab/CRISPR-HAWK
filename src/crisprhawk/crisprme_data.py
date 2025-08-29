""" """

from .crisprhawk_error import CrisprHawkPrepareDataError
from .exception_handlers import exception_handler
from .utils import IUPAC_ENCODER

from typing import Tuple, List, Union

import re
import os


def is_pam_right(header: str) -> bool:
    # if pam column occurs before, pam occurs upstream guide
    sgrna_idx = header.find("sgRNA_sequence")
    pam_idx = header.find("pam")
    assert sgrna_idx != pam_idx
    return True if pam_idx < sgrna_idx else False


def read_guides_report(report: str, debug: bool) -> Tuple[bool, List[List[str]]]:
    try:
        with open(report, mode="r") as infile:
            # assess if pam occurs upstream guides
            right = is_pam_right(infile.readline().strip())
            fields = [line.strip().split()[:6] for line in infile]  # first 6 cols
    except IOError as e:
        exception_handler(
            CrisprHawkPrepareDataError,
            f"Reading report {report} failed",
            os.EX_IOERR,
            debug,
            e,
        )
    return right, fields


def _replacer(match: re.Match) -> Union[str, None]:
    key = match.group(1)
    return IUPAC_ENCODER.get(key, match.group(0))


def solve_pam(pamclass: str) -> str:
    return re.sub(r"\[([^\]]+)\]", _replacer, pamclass)  # type: ignore


def create_pam_file(
    pamclass: str, right: bool, guidelen: int, outdir: str, debug: bool
) -> None:
    pam = solve_pam(pamclass)  # retrieve original pam with iupac
    try:  # write pam file
        with open(os.path.join(outdir, f"{pam}.txt"), mode="w") as outfile:
            gseq = "N" * guidelen  # add as many Ns as the guide length
            pamfull = (
                f"{pam}{gseq}\t{-len(pam)}" if right else f"{gseq}{pam}\t{len(pam)}"
            )
            outfile.write(f"{pamfull}\n")  # write PAM
    except IOError as e:
        exception_handler(
            CrisprHawkPrepareDataError,
            f"PAM file creation failed for PAM {pamclass}",
            os.EX_IOERR,
            debug,
            e,
        )
    if (
        not os.path.isfile(os.path.join(outdir, f"{pam}.txt"))
        or os.stat(os.path.join(outdir, f"{pam}.txt")).st_size <= 0
    ):
        exception_handler(
            CrisprHawkPrepareDataError,
            f"PAM file empty or not created for PAM {pamclass}",
            os.EX_IOERR,
            debug,
        )


def create_guide_files(
    guides_data: List[List[str]], right: bool, outdir: str, debug: bool
) -> None:
    pidx = 3 if right else 4  # assess pam column index
    gidx = 4 if right else 3  # assess guide sequence column index
    for gdata in guides_data:
        gfname = os.path.join(
            outdir, f"{gdata[0]}_{gdata[1]}_{gdata[2]}_{gdata[gidx]}_{gdata[pidx]}.txt"
        )
        try:
            with open(gfname, mode="w") as outfile:
                pamseq = "N" * len(gdata[pidx])  # add as many Ns as the pam length
                guidefull = (
                    f"{pamseq}{gdata[gidx]}" if right else f"{gdata[gidx]}{pamseq}"
                )
                outfile.write(f"{guidefull}\n")
        except IOError as e:
            exception_handler(
                CrisprHawkPrepareDataError,
                f"Guide file creation failed for guide {gfname}",
                os.EX_IOERR,
                debug,
                e,
            )
        if not os.path.isfile(gfname) or os.stat(gfname).st_size <= 0:
            exception_handler(
                CrisprHawkPrepareDataError,
                f"Guide file empty or not created for guide {gfname}",
                os.EX_IOERR,
                debug,
            )


def prepare_data_crisprme(
    report: str, create_pam: bool, outdir: str, debug: bool
) -> None:
    # read input guides report
    isright, guides_data = read_guides_report(report, debug)
    if create_pam:  # pam file creation requested
        gdata = guides_data[0]  # retrieve guide's data
        gidx = 4 if isright else 3  # assess guide sequence column index
        create_pam_file(gdata[5], isright, len(gdata[gidx]), outdir, debug)
    # create guide file for each guide in the report
    create_guide_files(guides_data, isright, outdir, debug)
