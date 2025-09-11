""" """

from .crisprhawk_error import (
    CrisprHawkOffTargetsError,
    CrisprHawkCfdScoreError,
    CrisprHawkElevationScoreError,
)
from .config_crispritz import CrispritzConfig, CRISPRITZ
from .exception_handlers import exception_handler
from .utils import (
    print_verbosity,
    suppress_stdout,
    suppress_stderr,
    VERBOSITYLVL,
)
from .region_constructor import PADDING
from .scores.cfdscore.cfdscore import load_mismatch_pam_scores
from .scores.crisprhawk_scores import elevation
from .offtarget import Offtarget
from .bedfile import BedAnnotation
from .region import Region
from .guide import Guide
from .pam import PAM, SPCAS9, XCAS9

from typing import List, Tuple, Dict
from time import time

import pandas as pd
import numpy as np

import subprocess
import os


# off-targets report column names
OTREPCNAMES = [
    "chrom",  # 0
    "position",  # 1
    "strand",  # 2
    "grna",  # 3
    "spacer",  # 4
    "pam",  # 5
    "mm",  # 6
    "bulge_size",  # 7
    "bulg_type",  # 8
    "cfd",  # 9
    "elevation",  # 10
]


def _write_guides_file(
    guides: List[Guide], pam: PAM, crispritz_dir: str, verbosity: int, debug: bool
) -> str:
    guides_fname = os.path.join(crispritz_dir, "guides.txt")  # guides file
    print_verbosity(
        f"Creating guides file for off-target estimation", verbosity, VERBOSITYLVL[3]
    )
    pamseq = "N" * len(pam)  # pam sequence in guide
    try:
        with open(guides_fname, mode="w") as outfile:
            for guide in guides:  # format each guide for crispritz input
                guide_f = (
                    f"{pamseq}{guide.guide.upper()}"
                    if guide.right
                    else f"{guide.guide}{pamseq}"
                )
                outfile.write(f"{guide_f}\n")  # write guide
    except (FileNotFoundError, IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing crispritz guide file {guides_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert os.stat(guides_fname).st_size > 0
    return guides_fname


def _write_pam_file(
    pam: PAM,
    guidelen: int,
    right: bool,
    crispritz_dir: str,
    verbosity: int,
    debug: bool,
) -> str:
    pam_fname = os.path.join(crispritz_dir, "pam.txt")  # pam file
    print_verbosity(
        f"Creating PAM file for off-targets estimation", verbosity, VERBOSITYLVL[3]
    )
    try:
        with open(pam_fname, mode="w") as outfile:
            gseq = "N" * guidelen  # guide sequence in pam
            pam_f = f"{pam}{gseq} {-len(pam)}" if right else f"{gseq}{pam} {len(pam)}"
            outfile.write(f"{pam_f}\n")  # write pam
    except (FileNotFoundError, IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing crispritz PAM file {pam_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert os.stat(pam_fname).st_size > 0
    return pam_fname


def _prepare_input_data(
    crispritz_config: CrispritzConfig,
    guides: List[Guide],
    pam: PAM,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> Tuple[str, str]:
    if crispritz_config.outdir == ".crispritz_targets":
        # create hidden folder within output directory
        crispritz_dir = os.path.join(outdir, crispritz_config.outdir)
    else:  # use directory reported in the config file
        crispritz_dir = crispritz_config.outdir
    if not os.path.isdir(crispritz_dir):  # stores crispritz targets
        os.makedirs(crispritz_dir)
    # create guides and pam files
    guides_fname = _write_guides_file(guides, pam, crispritz_dir, verbosity, debug)
    pam_fname = _write_pam_file(
        pam, guides[0].guidelen, guides[0].right, crispritz_dir, verbosity, debug
    )
    return guides_fname, pam_fname


def _format_targets_prefix(
    region: Region, pam: PAM, guidelen: int, crispritz_dir: str
) -> str:
    return os.path.join(
        crispritz_dir,
        f"{region.contig}_{region.start + PADDING}_{region.stop - PADDING}_{pam.pam}_{guidelen}",
    )


def search(
    crispritz_config: CrispritzConfig,
    crispritz_index: str,
    region: Region,
    pam: PAM,
    guidelen: int,
    guide_fname: str,
    pam_fname: str,
    mm: int,
    bdna: int,
    brna: int,
    threads: int,
    crispritz_dir: str,
    verbosity: int,
    debug: bool,
) -> str:
    start = time()
    print_verbosity("Searching off-targets with CRISPRitz", verbosity, VERBOSITYLVL[2])
    # define crispritz command
    targets_prefix = _format_targets_prefix(region, pam, guidelen, crispritz_dir)
    crispritz_run = (
        f"{crispritz_config.conda} run -n {crispritz_config.env_name} {CRISPRITZ} "
        f"search {crispritz_index} {pam_fname} {guide_fname} {targets_prefix} -mm "
        f"{mm} -bDNA {bdna} -bRNA {brna} -th {threads} -r"
    )
    try:  # run crispritz to search for off-targets
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(
                crispritz_run,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    except subprocess.CalledProcessError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to search offtargets for guides in {guide_fname}",
            os.EX_DATAERR,
            debug,
            e,
        )
    targets_fname = f"{targets_prefix}.targets.txt"
    assert os.stat(targets_fname).st_size > 0
    print_verbosity(
        f"Off-targets search with CRISPRitz completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return targets_fname


def _read_offtargets(
    crispritz_targets_file: str, pam: PAM, right: bool, debug: bool
) -> List[Offtarget]:
    offtargets = []  # off-targets list
    try:
        with open(crispritz_targets_file, mode="r") as infile:
            infile.readline()  # skip header
            for line in infile:  # read crispritz report
                offtargets.append(Offtarget(line, pam.pam, right, debug))  # parse line
    except (FileNotFoundError, IOError, Exception) as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed retrieving CRISPRitz off-targets in {crispritz_targets_file}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return offtargets


def _compute_cfd_score(
    offtargets: List[Offtarget], verbosity: int, debug: bool
) -> List[Offtarget]:
    print_verbosity(
        f"Computing CFD score for {len(offtargets)} off-targets",
        verbosity,
        VERBOSITYLVL[2],
    )
    start = time()
    mmscores, pamscores = load_mismatch_pam_scores(debug)
    try:
        for ot in offtargets:
            ot.compute_cfd(mmscores, pamscores)
    except (ValueError, Exception) as e:
        exception_handler(
            CrisprHawkCfdScoreError,
            "Off-targets CFD scoring failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"CFD score computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return offtargets


def _compute_elevation_score(
    offtargets: List[Offtarget], verbosity: int, debug: bool
) -> List[Offtarget]:
    print_verbosity(
        f"Computing Elevation score for {len(offtargets)} off-targets",
        verbosity,
        VERBOSITYLVL[2],
    )
    start = time()
    # create wildtype and offtarget lists
    wildtypes_list, offtargets_list = [], []
    offtargets_scored, offtargets_filt = [], []
    for ot in offtargets:  # if bulge present skip and keep NA on elevation
        if not ("-" in ot.grna or "-" in ot.spacer):
            wildtypes_list.append(ot.grna.upper())  # wildtypes
            offtargets_list.append(ot.spacer.upper())  # offtargets
            offtargets_scored.append(ot)  # offtargets scored with elevation
        else:
            offtargets_filt.append(ot)
    try:
        for i, score in enumerate(elevation(wildtypes_list, offtargets_list)):
            offtargets_scored[i].set_elevation(score)

    except (ValueError, Exception) as e:
        exception_handler(
            CrisprHawkElevationScoreError,
            "Off-targets Elevation scoring failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"Elevation score computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return offtargets_filt + offtargets_scored


def report_offtargets(
    crispritz_targets_file: str,
    region: Region,
    pam: PAM,
    guidelen: int,
    right: bool,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> List[Offtarget]:
    # write haplotypes table to file
    start = time()  # track haplotypes table writing time
    offtargets = _read_offtargets(crispritz_targets_file, pam, right, debug)
    if pam.cas_system in [SPCAS9, XCAS9]:  # compute CFD score
        offtargets = _compute_cfd_score(offtargets, verbosity, debug)
    if guidelen + len(pam) == 23 and not right:  # compute elevation score
        offtargets = _compute_elevation_score(offtargets, verbosity, debug)
    print_verbosity("Writing off-targets report", verbosity, VERBOSITYLVL[1])
    report_fname = os.path.join(
        outdir,
        f"offtargets_{region.contig}_{region.start + PADDING}_{region.stop - PADDING}.tsv",
    )
    try:
        with open(report_fname, mode="w") as outfile:
            outfile.write("\t".join(OTREPCNAMES) + "\n")
            outfile.write("\n".join([ot.report_line() for ot in offtargets]))
        ot_table = pd.read_csv(report_fname, sep="\t")
        ot_table = ot_table.sort_values(OTREPCNAMES[:2])
        ot_table.to_csv(report_fname, sep="\t", index=False, na_rep="NA")
    except OSError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed writing off-targets report for region {region}",
            os.EX_IOERR,
            debug,
            e,
        )
    print_verbosity(
        f"Off-targets report written in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[2],
    )
    return offtargets


def _calculate_offtargets_map(
    offtargets: List[Offtarget], guides: List[Guide]
) -> Dict[str, List[Offtarget]]:
    otmap = {g.guide.upper(): [] for g in guides}  # offtargets map
    for ot in offtargets:
        # add each spacer to the corresponding guide (no pam)
        otmap[ot.grna_.upper().replace("-", "")].append(ot)
    return otmap


def _calculate_global_cfd(offtargets: List[Offtarget], verbosity: int) -> float:
    print_verbosity("Computing guide global CFD", verbosity, VERBOSITYLVL[2])
    start = time()
    cfds = [0 if ot.cfd == "NA" else float(ot.cfd) for ot in offtargets]
    print_verbosity(
        f"Global CFD computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return 100 / (100 + sum(cfds))


def _calculate_global_elevation(offtargets: List[Offtarget], verbosity: int) -> float:
    print_verbosity("Computing guide global Elevation", verbosity, VERBOSITYLVL[2])
    start = time()
    elevations = [
        0 if ot.elevation == "NA" else float(ot.elevation) for ot in offtargets
    ]
    print_verbosity(
        f"Global Elevation computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return 100 / (100 + sum(elevations))


def annotate_guides_offtargets(
    offtargets: List[Offtarget], guides: List[Guide], verbosity: int
) -> List[Guide]:
    otmap = _calculate_offtargets_map(offtargets, guides)
    for guide in guides:
        guide.set_offtargets(len(otmap[guide.guide]))  # set off-targets number
        guide.set_cfd(
            _calculate_global_cfd(otmap[guide.guide], verbosity)
        )  # set global CFD
        guide.set_elevation(
            _calculate_global_elevation(otmap[guide.guide], verbosity)
        )  # set global elevation
    return guides


def search_offtargets(
    guides: List[Guide],
    pam: PAM,
    crispritz_index: str,
    region: Region,
    crispritz_config: CrispritzConfig,
    mm: int,
    bdna: int,
    brna: int,
    guidelen: int,
    right: bool,
    threads: int,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    print_verbosity(
        f"Estimating off-targets for {len(guides)} guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()
    guides_fname, pam_fname = _prepare_input_data(
        crispritz_config, guides, pam, outdir, verbosity, debug
    )
    # search offtargets with crispritz
    targets_fname = search(
        crispritz_config,
        crispritz_index,
        region,
        pam,
        guides[0].guidelen,
        guides_fname,
        pam_fname,
        mm,
        bdna,
        brna,
        threads,
        os.path.dirname(guides_fname),
        verbosity,
        debug,
    )
    offtargets = report_offtargets(
        targets_fname, region, pam, guidelen, right, outdir, verbosity, debug
    )
    guides = annotate_guides_offtargets(offtargets, guides, verbosity)
    print_verbosity(
        f"Off-targets estimation completed in {(time() - start):.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides
