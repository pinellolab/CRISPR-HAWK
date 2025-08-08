""" """

from .crisprhawk_error import CrisprHawkOffTargetsError
from .exception_handlers import exception_handler
from .utils import (
    print_verbosity,
    suppress_stdout,
    suppress_stderr,
    create_temp_folder,
    remove_folder,
    IUPACTABLE,
    VERBOSITYLVL,
    DNA,
)
from .scores.cfdscore.cfdscore import load_mismatch_pam_scores
from .offtarget import Offtarget
from .bedfile import BedAnnotation
from .region import Region
from .guide import Guide
from .pam import PAM

from dataclasses import dataclass
from typing import List, Tuple, NamedTuple
from pybedtools import BedTool
from time import time

import multiprocessing
import subprocess
import os


MAXMM = 4  # maximum number of mismatches allowed when searching offtargets
MAXENTRIES = 2000000  # maximum entries (bwa param)
MAXOCC = 60000  # maximum number of occurrences (bwa param)
OTREPCNAMES = [
    "guide",
    "chrom",
    "start",
    "stop",
    "strand",
    "spacer",
    "pam",
    "mm",
    "cfd",
    "functional_annotation",
    "gene_annotation",
]


@dataclass
class GuideProcessingArgs:
    guide: Guide
    pam: PAM
    genome: str
    functional_annotation: str
    gene_annotation: str
    write_offtargets_report: bool
    verbosity: int
    debug: bool


class GuideProcessingResult(NamedTuple):
    guide: Guide
    offtarget_report_lines: List[str]
    num_offtargets: int
    cfd_score: float


def create_guide_fa(guide: Guide, workingdir: str, debug: bool) -> str:
    fname = os.path.join(workingdir, f"{guide.guide_id}.fa")
    try:  # create dummy fasta containing guide's  sequence for alignment with bwa
        with open(fname, mode="w") as outfile:
            outfile.write(f">{guide.guide_id}\n{guide.guide}\n")
    except OSError as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            "Failed to create dummy guide FASTA file",
            os.EX_IOERR,
            debug,
            e,
        )
    return fname


def search(guide_fa: str, genome: str, guide_len: int, debug: bool) -> str:
    guide_sa = f"{guide_fa}.sa"  # sa guide file
    if not os.path.isfile(f"{genome}.sa"):
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Genome index {genome}.sa not found",
            os.EX_DATAERR,
            debug,
        )
    try:
        # run bwa to search for off-targets
        cmd = f"bwa aln -o 0 -m {MAXENTRIES} -k {MAXMM} -n {MAXMM} -N -l {guide_len} {genome} {guide_fa} > {guide_sa}"
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
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
    script_bin = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "scripts/xa2multi.pl"
    )
    try:
        cmd = f"bwa samse -n {MAXOCC} {genome} {guide_sa} {guide_fa} | perl {script_bin} > {guide_sam}"
        with suppress_stdout(), suppress_stderr():
            # suppress stdout and stderr to avoid cluttering the output
            subprocess.check_call(
                cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
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
            x1score = int(value)  # alternative alignment score
        elif tagname == "XA":
            hasalt = True  # alternative alignments are present
    assert mm is not None  # mismatches must be at least 0
    return (
        mm,
        x0score + x1score,
        hasalt,
    )  # return number of mismatches, total alignment score and alternative alignments flag


def adjust_target_positions(
    position: int, right: bool, strand: str, pamlen: int, guidelen: int
) -> Tuple[int, int]:
    startpos = position - 1  # convert to 0-based index
    endpos = startpos + guidelen  # calculate end position based on guide length
    if (right and strand == "+") or (not right and strand == "-"):
        startpos -= pamlen
    else:
        endpos += pamlen
    return startpos, endpos  # return adjusted start and end positions


def _read_sam_line(
    line: str, guidelen: int, pamlen: int, guide_pos: int, right: bool, debug: bool
) -> str:
    fields = line.strip().split()  # retrieve individual fields
    if fields[5] == "*":
        return ""  # skip unmapped (empty cigar)
    # retrieve mm, alignment score, alt alignments
    mm, alnscore, hasalt = parse_tags(fields[11:])
    isrep = "1" if hasalt and alnscore > 1 else "0"  # repeated alignement
    strand = "-" if (int(fields[1]) & 16) == 16 else "+"  # determine strand
    if fields[5] != f"{guidelen}M":
        exception_handler(
            CrisprHawkOffTargetsError,
            f"Unexpected CIGAR format {fields[5]} in line {line.strip()}",
            os.EX_DATAERR,
            debug,
        )
    startpos, endpos = adjust_target_positions(
        int(fields[3]), right, strand, pamlen, guidelen
    )  # adjust target's start and stop positions
    if startpos < 0 or startpos == (
        guide_pos - 1
    ):  # skip off-targets at sequence limits
        return ""
    # write to bed file: chrom, start, end, name, score, strand, mm, isrep
    targetname = "|".join(
        [
            fields[0],
            fields[2],
            str(startpos),
            str(endpos),
            strand,
            str(mm),
            str(alnscore),
            isrep,
        ]
    )
    return "\t".join(
        [fields[2], str(startpos), str(endpos), targetname, str(mm), strand]
    )


def generate_bed(
    guide_sam: str, guide_pos: int, pamlen: int, guidelen: int, right: bool, debug: bool
) -> str:
    guide_bed = f"{guide_sam}.bed"
    try:  # retrieve bed file for found off-targets from bwa sam
        with open(guide_sam, mode="r") as infile, open(guide_bed, mode="w") as outfile:
            for line in infile:
                if line.startswith("@"):  # skip comments (sam)
                    continue
                if not (
                    bedline := _read_sam_line(
                        line, guidelen, pamlen, guide_pos, right, debug
                    )
                ):
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


def sam2bed(
    guide_fa: str,
    guide_sa: str,
    genome: str,
    position: int,
    pamlen: int,
    guidelen: int,
    right: bool,
    debug: bool,
) -> str:
    # generate sam file from bwa alignment
    guide_sam = generate_samse(guide_fa, guide_sa, genome, debug=False)
    # generate bed file from sam
    guide_bed = generate_bed(guide_sam, position, pamlen, guidelen, right, debug)
    try:
        cmd = f"sort -k1,1 -k2,2n {guide_bed} | bedClip stdin {genome}.sizes stdout > {guide_bed}.sorted"
        subprocess.check_call(
            cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
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
    bed = BedTool(guide_bed)  # load bed file
    sequences = bed.sequence(
        fi=genome, s=True, name=True  # type: ignore
    )  # extract offtargets sequences 
    # construct a list of offtarget sequences
    sequences = open(sequences.seqfn).read().split()
    return [(sequences[i], sequences[i + 1]) for i in range(0, len(sequences), 2)]


def matchpam(pam: str, otpam: str) -> bool:
    assert len(pam) == len(otpam)
    for i, nt in enumerate(pam):
        if otpam[i] not in IUPACTABLE[nt]:
            return False
    return True


def filter_offtargets(
    sequences: List[Tuple[str, str]], pam: PAM, pamlen: int, right: bool, debug: bool
) -> List[Offtarget]:
    offtargets = []
    for otannotation, offtarget in sequences:
        # skip off-targets containing Ns
        if DNA[4] in offtarget.upper():
            continue
        # parse annotation fields: format depends on alignment tool output (bwa)
        annotation_fields = otannotation.split("|")
        # check if alignment should be skipped due to repetitive regions
        is_repetitive = bool(int(annotation_fields[7][:1]))
        if is_repetitive:
            continue  # skip repetitive alignments to avoid false positives
        # extract PAM sequence based on orientation
        otpam = offtarget[:pamlen] if right else offtarget[-pamlen:]
        otsequence = offtarget[pamlen:] if right else offtarget[:-pamlen]
        # validate PAM sequence against expected pattern
        if not matchpam(pam.pam, otpam.upper()):
            continue
        # convert 0-based BED coordinates to 1-based genomic coordinates
        start = int(annotation_fields[2]) + 1
        stop = int(annotation_fields[3]) + 1
        offtargets.append(
            Offtarget(
                annotation_fields[1],
                start,
                stop,
                annotation_fields[4],
                otsequence,
                otpam,
                int(annotation_fields[5]),
                debug,
            )
        )
    return offtargets


def retrieve_offtargets(
    guide_bed: str, genome: str, pam: PAM, pamlen: int, right: bool, debug: bool
) -> List[Offtarget]:
    offtargets = offtarget_sequence(guide_bed, genome)  # retrieve offtargets sequence
    return filter_offtargets(offtargets, pam, pamlen, right, debug)


def annotate_offtargets(
    wildtype: Guide,
    offtargets: List[Offtarget],
    functional_annotation: str,
    gene_annotation: str,
    verbosity: int,
    debug: bool,
) -> List[Offtarget]:
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
            offtarget.annotate_functional(funcann)  # type: ignore
        if annotate_gene:  # annotate off-target gene function
            offtarget.annotate_gene(geneann)  # type: ignore
    return offtargets


def compute_cfd_round(offtargets: List[Offtarget]) -> float:
    cfdscores = [ot.cfd for ot in offtargets]
    cfd = 100 / (100 + sum(cfdscores))
    return int(round(cfd * 100))


# TODO: compute elevation score


def report_offtargets(
    offtargets: List[str], region: Region, outdir: str, verbosity: int, debug: bool
) -> None:
    # write haplotypes table to file
    print_verbosity("Writing off-targets report", verbosity, VERBOSITYLVL[1])
    start = time()  # track haplotypes table writing time
    report_fname = os.path.join(
        outdir, f"offtargets_{region.contig}_{region.start}_{region.stop}.tsv"
    )
    try:
        with open(report_fname, mode="w") as outfile:
            outfile.write("\t".join(OTREPCNAMES) + "\n")
            outfile.write("\n".join(offtargets))
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


def _process_guide_offtargets(
    args: GuideProcessingArgs,
) -> Tuple[List[Offtarget], int, float]:
    workingdir = create_temp_folder(args.guide.guide_id)
    try:
        print_verbosity(
            f"Searching off-targets for guide {args.guide.guidepam}",
            args.verbosity,
            VERBOSITYLVL[3],
        )
        guide_fa = create_guide_fa(args.guide, workingdir, args.debug)  # guide fasta
        guide_sa = search(
            guide_fa, args.genome, args.guide.guidelen, args.debug
        )  # search offtargets
        guide_bed = sam2bed(
            guide_fa,
            guide_sa,
            args.genome,
            args.guide.start,
            args.guide.pamlen,
            args.guide.guidelen,
            args.guide.right,
            args.debug,
        )  # convert sam to bed
        offtargets_guide = retrieve_offtargets(
            guide_bed,
            args.genome,
            args.pam,
            args.guide.pamlen,
            args.guide.right,
            args.debug,
        )
        offtargets_guide = annotate_offtargets(
            args.guide,
            offtargets_guide,
            args.functional_annotation,
            args.gene_annotation,
            args.verbosity,
            args.debug,
        )
        # compute cfd and offtargets number for each guide
        num_offtargets = len(offtargets_guide) if offtargets_guide else -1
        cfd_score = compute_cfd_round(offtargets_guide) if offtargets_guide else 0.0
        return offtargets_guide, num_offtargets, cfd_score
    finally:  # remove temporary working directory
        remove_folder(workingdir)


def process_single_guide(args: GuideProcessingArgs) -> GuideProcessingResult:
    try:
        offtargets_guide, num_offtargets, cfd_score = _process_guide_offtargets(args)
        # generate report lines if needed
        offtargets_report_lines = []
        if args.write_offtargets_report and offtargets_guide:
            offtargets_report_lines = [
                offtarget.report_line(args.guide) for offtarget in offtargets_guide
            ]
        return GuideProcessingResult(
            guide=args.guide,
            offtarget_report_lines=offtargets_report_lines,
            num_offtargets=num_offtargets,
            cfd_score=cfd_score,
        )
    except Exception as e:
        raise exception_handler(
            CrisprHawkOffTargetsError,
            f"Off-target search failed for guide {args.guide.guidepam}",
            os.EX_DATAERR,
            args.debug,
            e,
        )


def _prepare_processing_args(
    guides: List[Guide],
    pam: PAM,
    genome: str,
    functional_annotation: str,
    gene_annotation: str,
    write_offtargets_report: bool,
    verbosity: int,
    debug: bool,
) -> List[GuideProcessingArgs]:
    return [
        GuideProcessingArgs(
            guide=guide,
            pam=pam,
            genome=genome,
            functional_annotation=functional_annotation,
            gene_annotation=gene_annotation,
            write_offtargets_report=write_offtargets_report,
            verbosity=verbosity,
            debug=debug,
        )
        for guide in guides
    ]


def _execute_parallel_processing(
    argslist: List[GuideProcessingArgs], threads: int
) -> List[GuideProcessingResult]:
    try:
        with multiprocessing.Pool(processes=threads) as pool:
            async_result = pool.map_async(process_single_guide, argslist)
            return async_result.get()
    except Exception as e:
        exception_handler(
            CrisprHawkOffTargetsError,
            "Off-target search failed",
            os.EX_DATAERR,
            argslist[0].debug,
            e,
        )


def _update_guides_from_results(
    results: List[GuideProcessingResult],
) -> Tuple[List[Guide], List[str]]:
    updated_guides, offtargets_lines = [], []
    for res in results:  # update guide with offtarget info and collect offtarget lines
        res.guide.set_offtargets(res.num_offtargets)
        res.guide.set_cfd(res.cfd_score)
        updated_guides.append(res.guide)
        offtargets_lines += res.offtarget_report_lines
    return updated_guides, offtargets_lines


def search_offtargets(
    guides: List[Guide],
    pam: PAM,
    genome: str,
    region: Region,
    functional_annotation: str,
    gene_annotation: str,
    write_offtargets_report: bool,
    threads: int,
    outdir: str,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    print_verbosity(
        f"Estimating off-targets for {len(guides)} guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()
    # prepare arguments for parallel processing
    argslist = _prepare_processing_args(
        guides,
        pam,
        genome,
        functional_annotation,
        gene_annotation,
        write_offtargets_report,
        verbosity,
        debug,
    )
    # execute parallel processing
    results = _execute_parallel_processing(argslist, threads)
    # update guides and collect report data
    updated_guides, offtarget_lines = _update_guides_from_results(results)
    # write off-targets report if requested
    if write_offtargets_report and offtarget_lines:
        report_offtargets(offtarget_lines, region, outdir, verbosity, debug)
    print_verbosity(
        f"Off-targets estimation completed in {(time() - start):.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return updated_guides
