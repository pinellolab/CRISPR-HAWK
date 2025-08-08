""" """

from .crisprhawk_error import CrisprHawkOffTargetsError
from .exception_handlers import exception_handler
from .utils import print_verbosity, suppress_stdout, suppress_stderr, IUPACTABLE, VERBOSITYLVL
from .scores.cfdscore.cfdscore import load_mismatch_pam_scores
from .offtarget import Offtarget
from .bedfile import BedAnnotation
from .region import Region
from .guide import Guide
from .pam import PAM

from typing import List, Tuple, Dict
from pybedtools import BedTool
from time import time

import numpy as np

import tempfile
import subprocess
import os

from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import functools

MAXMM = 4  # maximum number of mismatches allowed when searching offtargets
MAXENTRIES = 2000000  # maximum entries (bwa param)
MAXOCC = 60000  # maximum number of occurrences (bwa param)
OTREPCNAMES = ["guide", "chrom", "start", "stop", "strand", "spacer", "pam", "mm", "cfd", "functional_annotation", "gene_annotation"]

def create_guide_fa(guide: Guide, workingdir: str, debug: bool) -> str:
    fname = os.path.join(workingdir, f"{guide.guide_id}.fa")
    try:  # create dummy fasta containing guide's  sequence for alignment with bwa
        with open(fname, mode="w") as outfile:
            outfile.write(f">{guide.guide_id}\n{guide.guide}\n")
    except OSError as e:
        exception_handler(CrisprHawkOffTargetsError, "Failed to create dummy guide FASTA file", os.EX_IOERR, debug, e)
    return fname

def search(guide_fa: str, genome: str, guide_len: int, debug: bool) -> str:
    guide_sa = f"{guide_fa}.sa"  # sa guide file
    if not os.path.isfile(f"{genome}.sa"):
        exception_handler(CrisprHawkOffTargetsError, f"Genome index {genome}.sa not found", os.EX_DATAERR, debug)
    try:
        # run bwa to search for off-targets
        cmd = f"bwa aln -o 0 -m {MAXENTRIES} -k {MAXMM} -n {MAXMM} -N -l {guide_len} {genome} {guide_fa} > {guide_sa}"
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to search offtargets for guide {guide_fa}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return guide_sa  # return sa file with off-targets

def generate_samse(guide_fa: str, guide_sa: str, genome: str, debug: bool) -> str:
    guide_sam = f"{guide_fa}.sam"  # sam guide file
    # perl script to construct feasible output
    script_bin = os.path.join(os.path.abspath(os.path.dirname(__file__)), "scripts/xa2multi.pl")
    try:
        cmd = f"bwa samse -n {MAXOCC} {genome} {guide_sa} {guide_fa} | perl {script_bin} > {guide_sam}"
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to run bwa samse for guide {guide_fa}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return guide_sam  # return sam file with off-targets

def parse_tags(tags: List[str]) -> Tuple[int, int, bool]:
    x0score, x1score = 0, 0  # initialize alignment scores
    mm = None  # initialize number of mismatches
    hasalt = False  # flag to check if alternative alignments are present
    for tag in tags:
        tagname, dtype, value = tag.split(":")  # retrieve tag name, type and value
        if tagname == "NM":
            mm = int(value)  # number of mismatches
        elif tagname == "X0":
            x0score = int(value)  # primary alignment score
        elif tagname == "X1":
            x1score = int(value)   # alternative alignment score
        elif tagname == "XA":
            hasalt = True  # alternative alignments are present
    assert mm is not None  # mismatches must be at least 0
    return mm, x0score + x1score, hasalt  # return number of mismatches, total alignment score and alternative alignments flag

def adjust_target_positions(position: int, right: bool, strand: str, pamlen: int, guidelen: int) -> Tuple[int, int]:
    startpos = position - 1  # convert to 0-based index
    endpos = startpos + guidelen  # calculate end position based on guide length
    if (right and strand == "+") or (not right and strand == "-"):
        startpos -= pamlen
    else:
        endpos += pamlen
    return startpos, endpos  # return adjusted start and end positions

def _read_sam_line(line: str, guidelen: int, pamlen: int, guide_pos: int, right: bool, debug: bool) -> str:
    fields = line.strip().split()  # retrieve individual fields
    if fields[5] == "*":
        return ""  # skip unmapped (empty cigar)
    # retrieve mm, alignment score, alt alignments
    mm, alnscore, hasalt = parse_tags(fields[11:])
    isrep = "1" if hasalt and alnscore > 1 else "0"  # repeated alignement
    strand = "-" if (int(fields[1]) & 16) == 16 else "+"  # determine strand
    if fields[5] != f"{guidelen}M":
        exception_handler(CrisprHawkOffTargetsError, f"Unexpected CIGAR format {fields[5]} in line {line.strip()}", os.EX_DATAERR, debug)
    startpos, endpos = adjust_target_positions(int(fields[3]), right, strand, pamlen, guidelen)  # adjust target's start and stop positions
    if startpos < 0 or startpos == (guide_pos - 1):  # skip off-targets at sequence limits
        return ""
    # write to bed file: chrom, start, end, name, score, strand, mm, isrep
    targetname = "|".join([fields[0], fields[2], str(startpos), str(endpos), strand, str(mm), str(alnscore), isrep])
    return "\t".join([fields[2], str(startpos), str(endpos), targetname, str(mm), strand])


def generate_bed(guide_sam: str, guide_pos: int, pamlen: int, guidelen: int, right: bool, debug: bool) -> str:
    guide_bed = f"{guide_sam}.bed"
    try:  # retrieve bed file for found off-targets from bwa sam
        with open(guide_sam, mode="r") as infile, open(guide_bed, mode="w") as outfile:
            for line in infile:
                if line.startswith("@"):  # skip comments (sam)
                    continue
                if not (bedline := _read_sam_line(line, guidelen, pamlen, guide_pos, right, debug)):
                    continue  # skip line if unmapped or mismapped
                outfile.write(f"{bedline}\n")
    except OSError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to convert SAM to BED for guide {guide_sam}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return guide_bed  # return bed file with off-targets

def sam2bed(guide_fa: str, guide_sa: str, genome: str, position: int, pamlen: int, guidelen: int, right: bool, debug: bool) -> str:
    # generate sam file from bwa alignment
    guide_sam = generate_samse(guide_fa, guide_sa, genome, debug=False)
    # generate bed file from sam
    guide_bed = generate_bed(guide_sam, position, pamlen, guidelen, right, debug)
    try:
        cmd = f"sort -k1,1 -k2,2n {guide_bed} | bedClip stdin {genome}.sizes stdout > {guide_bed}.sorted"
        subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Failed to sort and clip BED file {guide_bed}",
            os.EX_DATAERR,
            debug,
            e,
        )
    return f"{guide_bed}.sorted"

def offtarget_sequence(guide_bed: str, genome: str) -> List[Tuple[str, str]]:
    bed = BedTool(guide_bed) # load bed file
    sequences = bed.sequence(fi=genome, s=True, name=True)  # extract offtargets sequences # type: ignore
    # construct a list of offtarget sequences
    sequences = open(sequences.seqfn).read().split()
    return [(sequences[i], sequences[i + 1]) for i in range(0, len(sequences), 2)]

def matchpam(pam: str, otpam: str) -> bool:
    assert len(pam) == len(otpam)
    for i, nt in enumerate(pam):
        if otpam[i] not in IUPACTABLE[nt]:
            return False
    return True

def filter_offtargets(sequences: List[Tuple[str, str]], pam: PAM, pamlen: int, right: bool, debug: bool) -> List[Offtarget]:
    offtargets = []
    for otannotation, offtarget in sequences:
        # parse annotation fields: format depends on alignment tool output (bwa)
        annotation_fields = otannotation.split("|")
        # check if alignment should be skipped due to repetitive regions
        is_repetitive = bool(int(annotation_fields[7][:1]))
        if is_repetitive:
            continue  # skip repetitive alignments to avoid false positives
        # extract PAM sequence based on orientation
        otpam = (offtarget[:pamlen] if right else offtarget[-pamlen:])
        otsequence = (offtarget[pamlen:] if right else offtarget[:-pamlen])
        # validate PAM sequence against expected pattern
        if not matchpam(pam.pam, otpam.upper()):
            continue
        # convert 0-based BED coordinates to 1-based genomic coordinates
        start = int(annotation_fields[2]) + 1
        stop = int(annotation_fields[3]) + 1
        offtargets.append(Offtarget(annotation_fields[1], start, stop, annotation_fields[4], otsequence, otpam, int(annotation_fields[5]), debug))
    return offtargets


def retrieve_offtargets(guide_bed: str, genome: str, pam: PAM, pamlen: int, right: bool, debug: bool) -> List[Offtarget]:
    offtargets = offtarget_sequence(guide_bed, genome)  # retrieve offtargets sequence
    return filter_offtargets(offtargets, pam, pamlen, right, debug)

def annotate_offtargets(wildtype: Guide, offtargets: List[Offtarget], functional_annotation: str, gene_annotation: str, verbosity: int, debug: bool) -> List[Offtarget]:
    # load mismatch and pam scoring models (CFD score)
    mmscores, pamscores = load_mismatch_pam_scores(debug)
    annotate_functional = bool(functional_annotation)
    annotate_gene = bool(gene_annotation)
    if annotate_functional:
        funcann = BedAnnotation(functional_annotation, verbosity, debug)
    if annotate_gene:
        geneann = BedAnnotation(gene_annotation, verbosity, debug)
    for offtarget in offtargets:
        offtarget.compute_cfd(wildtype, mmscores, pamscores)
        if annotate_functional:  # annotate off-target functionally
            offtarget.annotate_functional(funcann) # type: ignore
        if annotate_gene:  # annotate off-target gene function
            offtarget.annotate_gene(geneann) # type: ignore
    return offtargets

def compute_cfd_round(offtargets: List[Offtarget]) -> float:
    cfdscores = [ot.cfd for ot in offtargets]
    cfd = 100 / (100 + sum(cfdscores))
    return int(round(cfd * 100))

# TODO: compute elevation score

def report_offtargets(offtargets: List[str], region: Region, outdir: str, verbosity: int, debug: bool) -> None:
    # write haplotypes table to file
    print_verbosity("Writing off-targets report", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes table writing time
    report_fname = os.path.join(outdir, f"offtargets_{region.contig}_{region.start}_{region.stop}.tsv")
    try:
        with open(report_fname, mode="w") as outfile:
            outfile.write("\t".join(OTREPCNAMES) + "\n")
            outfile.write("\n".join(offtargets))
    except OSError as e:
            exception_handler(CrisprHawkOffTargetsError, f"Failed writing off-targets report for region {region}", os.EX_IOERR, debug, e)  
    print_verbosity(
        f"Off-targets report written in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )


def search_offtargets(guides: List[Guide], pam: PAM, genome: str, region: Region, functional_annotation: str, gene_annotation: str, write_offtargets_report: bool, outdir: str, verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity(f"Estimating off-targets for {len(guides)} guides", verbosity, VERBOSITYLVL[3])
    start = time()
    offtargets = []  # list storing the all off-targets for writing ots report
    workingdir = tempfile.mkdtemp(prefix="crisprhawk_")  # create temp dir
    try:
        for guide in guides:
            print_verbosity(f"Searching off-targets for guide {guide.guidepam}", verbosity, VERBOSITYLVL[3])
            guide_fa = create_guide_fa(guide, workingdir, debug)  # guide fasta
            guide_sa = search(guide_fa, genome, guide.guidelen, debug)  # search offtargets
            guide_bed = sam2bed(guide_fa, guide_sa, genome, guide.start, guide.pamlen, guide.guidelen, guide.right, debug)  # convert sam to bed
            offtargets_guide = retrieve_offtargets(guide_bed, genome, pam, guide.pamlen, guide.right, debug)
            offtargets_guide = annotate_offtargets(guide, offtargets_guide, functional_annotation, gene_annotation, verbosity, debug)
            if write_offtargets_report:
                offtargets += [offtarget.report_line(guide) for offtarget in offtargets_guide]
            # set number of offtragets and guide's cfd
            if len(offtargets_guide) > 0:
                guide.set_offtargets(len(offtargets_guide))
                guide.set_cfd(compute_cfd_round(offtargets_guide))
            else:
                guide.set_offtargets(-1)  # type: ignore
                guide.set_cfd(0.)
            # clean folder from tmp files 
            tmpfiles = os.path.join(workingdir, f"{guide.guide_id}.*")
            subprocess.call(f"rm -r {tmpfiles}", shell=True)
    except Exception:
        subprocess.call(f"rm -r {workingdir}", shell=True)  # remove temp dir
        raise Exception("Off-target search failed")
    if write_offtargets_report:  # write offtargets report analyzed guides
        report_offtargets(offtargets, region, outdir, verbosity, debug)
    print_verbosity(
        f"Off-targets estimation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides



# TODO: off-targets search in parallel

# def _process_single_guide(guide: Guide, pam: PAM, genome: str, functional_annotation: str, gene_annotation: str, working_directory: str, verbosity: int, debug: bool) -> List[Offtarget]:
#     print_verbosity(f"Searching off-targets for guide {guide.guidepam}", verbosity, VERBOSITYLVL[3])
#     # create guide-specific FASTA file
#     guide_fasta = create_guide_fa(guide, working_directory, debug)
#     # search for off-targets using BWA
#     guide_alignment = search(guide_fasta, genome, guide.guidelen, debug)
#     # convert alignment results to BED format
#     guide_bed = sam2bed(guide_fasta, guide_alignment, genome, guide.start, guide.pamlen, guide.guidelen, guide.right, debug)
#     # clean folder from tmp files 
#     tmpfiles = os.path.join(working_directory, f"{guide.guide_id}.*")
#     subprocess.call(f"rm -r {tmpfiles}", shell=True)
#     # retrieve and filter off-target sequences
#     guide_offtargets = retrieve_offtargets(guide_bed, genome, pam, guide.pamlen, guide.right, debug)
#     # annotate off-targets with functional and gene information
#     annotated_offtargets = annotate_offtargets(guide, guide_offtargets, functional_annotation, gene_annotation, verbosity, debug)
#     return annotated_offtargets


# def _search_offtargets(guides: List[Guide], pam: PAM, genome: str, functional_annotation: str, gene_annotation: str, threads: int, verbosity: int, debug: bool) -> Dict[Guide, List[Offtarget]]:
#     print_verbosity(f"Estimating off-targets for {len(guides)} guides", verbosity, VERBOSITYLVL[3])
#     start_time = time()
#     threads = min(len(guides), cpu_count())  # determine optimal number of threads
#     # create temporary working directory
#     working_directory = tempfile.mkdtemp(prefix="crisprhawk_")
#     guide_offtargets = {}
#     try:  # create partial function with fixed parameters for parallel execution
#         process_guide = functools.partial(
#             _process_single_guide,
#             pam=pam,
#             genome=genome,
#             functional_annotation=functional_annotation,
#             gene_annotation=gene_annotation,
#             working_directory=working_directory,
#             verbosity=verbosity,
#             debug=debug
#         )
#         # execute guide processing in parallel
#         with ProcessPoolExecutor(max_workers=threads) as executor:
#             future_to_guide = {executor.submit(process_guide, guide): guide for guide in guides}
#             for future in as_completed(future_to_guide):  # collect results 
#                 guide = future_to_guide[future]
#                 try:
#                     offtargets = future.result()
#                     guide_offtargets[guide] = offtargets
#                     print_verbosity(f"Completed off-target search for guide {guide.guidepam}: {len(offtargets)} off-targets found", verbosity, VERBOSITYLVL[3])
#                 except Exception as e:
#                     raise        
#     except Exception as e:  # clean up temporary directory on error
#         subprocess.call(f"rm -r {working_directory}", shell=True)
#         exception_handler(
#             CrisprHawkOffTargetsError,
#             f"Off-target search failed: {str(e)}",
#             os.EX_DATAERR,
#             debug,
#             e
#         )
#     else:  # clean up temporary directory on success
#         subprocess.call(f"rm -r {working_directory}", shell=True)
#     print_verbosity(
#         f"Off-targets estimation completed in {time() - start_time:.2f}s",
#         verbosity,
#         VERBOSITYLVL[3]
#     )
#     return guide_offtargets


# # Alternative version that maintains the original function signature (no return value)
# def search_offtargets(
#     guides: List[Guide], 
#     pam: PAM, 
#     genome: str, 
#     functional_annotation: str, 
#     gene_annotation: str, 
#     verbosity: int, 
#     debug: bool,
# ) -> None:
#     """Original function signature version - processes guides and modifies them in-place."""
#     guide_offtargets = _search_offtargets(guides, pam, genome, functional_annotation, gene_annotation, 16, verbosity, debug)
#     # for guide, offtargets in guide_offtargets.items():
#     #     # Assuming Guide class has a method to store off-targets
#     #     if hasattr(guide, 'set_offtargets'):
#     #         guide.set_offtargets(offtargets)
