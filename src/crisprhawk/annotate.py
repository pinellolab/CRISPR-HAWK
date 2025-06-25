""" """

from .crisprhawk_error import CrisprHawkCfdScoreError, CrisprHawkAzimuthScoreError, CrisprHawkRs3ScoreError, CrisprHawkAnnotationError, CrisprHawkOffTargetsError
from .exception_handlers import exception_handler
from .scores import azimuth, cfdon, rs3
from .bedfile import BedAnnotation
from .guide import Guide, GUIDESEQPAD
from .utils import print_verbosity, flatten_list, suppress_stderr, suppress_stdout, VERBOSITYLVL, IUPACTABLE
from .region import Region
from .pam import PAM

from pybedtools import BedTool
from collections import defaultdict
from typing import List, Dict, Union, Tuple
from time import time

import numpy as np

import subprocess
import tempfile
import os

ANNDIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), "annotations")  # annotation data directory
MAXMM = 4  # maximum number of mismatches allowed when searching offtargets
MAXENTRIES = 2000000  # maximu entries (bwa param)
MAXOCC = 60000  # maximum number of occurrences (bwa param)


def reverse_guides(guides: List[Guide], verbosity: int) -> List[Guide]:
    # compute reverse complement sequence for guides occurring on reverse strand
    print_verbosity(
        "Reversing guides occurring on reverse strand", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # reversal start time
    for guide in guides:
        if guide.strand == 1:  # guide on reverse strand
            guide.reverse_complement()
    print_verbosity(
        f"Guides reversed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def azimuth_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    # create guides np.ndarray required by azimuth; each guide must have 4 nts
    # upstream the guide sequence, and 3 nts downstream the pam
    print_verbosity("Computing Azimuth score", verbosity, VERBOSITYLVL[3])
    start = time()  # azimuth score start time
    guides_seqs = np.array(
        [
            guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
            for guide in guides
        ]
    )
    try:  # compute azimuth scores
        azimuth_scores = azimuth(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkAzimuthScoreError, "Azimuth score calculation failed", os.EX_DATAERR, debug, e
        )
    assert len(azimuth_scores) == len(guides)  # should match
    for i, score in enumerate(azimuth_scores):
        guides[i].set_azimuth_score(score)  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def rs3_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity("Computing RS3 score", verbosity, VERBOSITYLVL[3])
    start = time()  # rs3 score start time
    guides_seqs = [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]
    try:  # compute azimuth scores
        rs3_scores = rs3(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkRs3ScoreError, "RS3 score calculation failed", os.EX_DATAERR, debug, e
        )
    assert len(rs3_scores) == len(guides)  # should match
    for i, score in enumerate(rs3_scores):
        guides[i].set_rs3_score(score)  # assign score to each guide
    print_verbosity(
        f"RS3 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def group_guides_position(
    guides: List[Guide], debug: bool
) -> Dict[str, Dict[int, Union[Guide, List[Guide]]]]:
    # dictionary to map guides to positions (0 -> ref; 1 -> alt)
    pos_guide = defaultdict(lambda: {0: None, 1: []})
    for guide in guides:
        poskey = f"{guide.start}_{guide.strand}"
        if guide.samples == "REF":  # reference guide
            if pos_guide[poskey][0] is not None:
                exception_handler(
                    CrisprHawkCfdScoreError,
                    f"Duplicate REF guide at position {guide.start}? CFDon calculation failed",
                    os.EX_DATAERR,
                    debug,
                )
            pos_guide[poskey][0] = guide  # type: ignore
        pos_guide[poskey][1].append(guide)  # type: ignore
    return pos_guide  # type: ignore


def cfdon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity("Computing CFDon score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    for _, gg in guide_groups.items():
        try:
            cfdon_scores = cfdon(gg[0], gg[1], debug)  # type: ignore
        except Exception as e:
            exception_handler(
                CrisprHawkCfdScoreError,
                "CFDon score calculation failed",
                os.EX_DATAERR,
                debug,
                e,
            )
        for i, score in enumerate(cfdon_scores):
            gg[1][i].set_cfdon_score(score)  # type: ignore
    # revert grouped guides by position into list
    guides = flatten_list([gg[1] for _, gg in guide_groups.items()])  # type: ignore
    print_verbosity(
        f"CFDon scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def annotate_variants(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    guides_lst = []  # reported guides
    print_verbosity(
        "Annotating variants occurring in guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # position calculation start time
    for guide in guides:
        guide_vars = set()  # variants occurring in variant
        for variant in guide.variants.split(","):
            if variant == "NA":  # no variants
                if guide_vars:
                    exception_handler(
                        ValueError, "Forbidden NA variant", os.EX_DATAERR, debug
                    )
                guide_vars.add("NA")
                break
            try:  # retrieve each variant position
                variant_position = int(variant.split("-")[1])
            except TypeError as e:
                exception_handler(
                    TypeError,
                    f"Variant {variant} seems to have a non int position",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
            # assess whether the snp occurs within the guide or is part of the haplotype
            if guide.start <= variant_position < guide.stop:
                guide_vars.add(variant)
        if guide_vars:
            guide.set_variants(",".join(sorted(guide_vars)))
            guides_lst.append(guide)
    print_verbosity(
        f"Variants annotated in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides_lst

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

def _samse(guide_fa: str, guide_sa: str, genome: str, debug: bool) -> str:
    guide_sam = f"{guide_fa}.sam"  # sam guide file
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


def _parse_tags(tags: List[str]) -> Tuple[int, int, bool]:
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


def _sam2bed(guide_sam: str, guide_pos: int, pamlen: int, guidelen: int, right: bool, debug: bool) -> str:
    guide_bed = f"{guide_sam}.bed"
    try:
        with open(guide_sam, mode="r") as infile, open(guide_bed, mode="w") as outfile:
            for line in infile:
                if line.startswith("@"):  # skip comments on sam
                    continue
                fields = line.strip().split()
                mm, alnscore, hasalt = _parse_tags(fields[11:])  # parse tags
                isrep = "1" if hasalt and alnscore > 1 else "0"  # repetitive alignment
                strand = "-" if (int(fields[1]) & 16) == 16 else "+"  # determine strand
                if fields[5] == "*":  # skip unmapped reads (empty cigar)
                    continue
                if fields[5] != f"{guidelen}M":
                    exception_handler(CrisprHawkOffTargetsError, f"Unexpected CIGAR format {fields[5]} in line {line.strip()}", os.EX_DATAERR, debug)
                startpos, endpos = adjust_target_positions(int(fields[3]), right, strand, pamlen, guidelen)  # adjust target positions
                if startpos < 0 or startpos == (guide_pos - 1):  # skip off-targets at sequence limits
                    continue
                # write to bed file: chrom, start, end, name, score, strand, mm, isrep
                targetname = "|".join([fields[0], fields[2], str(startpos), str(endpos), strand, str(mm), str(alnscore), isrep])
                outfile.write("\t".join([fields[2], str(startpos), str(endpos), targetname, str(mm), strand]) + "\n")
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
    guide_sam = _samse(guide_fa, guide_sa, genome, debug=False)  # get sam file
    guide_bed = _sam2bed(guide_sam, position, pamlen, guidelen, right, debug)
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
    sequences = bed.sequence(fi=genome, s=True, name=True)  # extract offtargets sequences
    # construct a list of offtarget sequences
    sequences = open(sequences.seqfn).read().split()
    return [(sequences[i], sequences[i + 1]) for i in range(0, len(sequences), 2)]


def matchpam(pam: str, otpam: str) -> bool:
    assert len(pam) == len(otpam)
    for i, nt in enumerate(pam):
        if otpam[i] not in IUPACTABLE[nt]:
            return False
    return True



def filter_offtargets(sequences: List[Tuple[str, str]], pam: PAM, pamlen: int, right: bool) -> List[Tuple[str, str, str, str, str]]:
    offtargets = []
    for otname, ot in sequences:
        fields = otname.split("|")  # retrieve offtarget annotation fields
        skip = True
        if bool(int(fields[7][:1])):  # ignore repetitive alignments
            skip = False
        else:
            pamot = ot[:pamlen] if right else ot[-pamlen:]
            skip = not matchpam(pam.pam, pamot.upper()) # check pam validity
        if not skip:
            offtargets.append((fields[1], fields[2], fields[3], fields[4], ot, fields[5]))
    return offtargets






def retrieve_offtargets(guide_bed: str, genome: str, pam: PAM, pamlen: int, right: bool) -> List[Tuple[str, str, str, str, str]]:
    offtargets = offtarget_sequence(guide_bed, genome)  # retrieve offtargets sequence
    return filter_offtargets(offtargets, pam, pamlen, right)



def search_offtargets(guides: List[Guide], pam: PAM, genome: str, debug: bool):
    workingdir = tempfile.mkdtemp(prefix="crisprhawk_")  # create temp dir
    print(workingdir)
    try:
        for guide in guides:
            guide_fa = create_guide_fa(guide, workingdir, debug)  # guide fasta
            guide_sa = search(guide_fa, genome, guide.guidelen, debug)
            guide_bed = sam2bed(guide_fa, guide_sa, genome, guide.start, guide.pamlen, guide.guidelen, guide.right, debug)  # convert sam to bed
            x = retrieve_offtargets(guide_bed, genome, pam, guide.pamlen, guide.right)
            print(len(x))
            print("\n".join(["\t".join(e) for e in x]))
            break

    except Exception:
        subprocess.call(f"rm -r {workingdir}", shell=True)  # remove temp dir
        raise Exception("Off-target search failed")


def _funcann(guide: Guide, bedannotation: BedAnnotation, contig: str, atype: str, idx: int) -> Guide:
    # fetch annotation features overlapped by input guide
    if not (annotation := bedannotation.fetch_features(contig, guide.start, guide.stop, idx)):
        annotation = "NA"  # no annotation feature overlapped by the input guide
    if atype == "gene":  # set gene annotation
        guide.set_gene_ann(annotation)
    else:  # set functional annotation
        guide.set_func_ann(annotation)
    return guide


def funcann_guides(guides: List[Guide], contig: str, annotation: str, atype: str, verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity("Starting guides functional annotation", verbosity, VERBOSITYLVL[3])
    start = time()  # functional annotation start time
    assert atype in {"func", "gene"}  # used to set the proper field in guides
    idx = 22 if atype == "gene" else 9
    bedann = BedAnnotation(annotation, verbosity, debug)  # load annotation bed
    try:
        guides_ann = [_funcann(guide, bedann, contig, atype, idx) for guide in guides]
    except Exception as e:
        exception_handler(CrisprHawkAnnotationError, "Guides functional annotation failed", os.EX_DATAERR, debug, e)
    assert len(guides) == len(guides_ann)
    print_verbosity(f"Guides functional annotation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3])
    return guides_ann


def annotate_guides(
    guides: Dict[Region, List[Guide]], functional_annotation: str, gene_annotation: str, pam: PAM, genome: str, estimate_offtargets: bool, verbosity: int, debug: bool
) -> Dict[Region, List[Guide]]:
    # annotate guides with scores, variants and adjust positions
    print_verbosity("Annotating guides", verbosity, VERBOSITYLVL[1])
    start = time()  # annotation start time
    for region, guides_list in guides.items():
        # compute reverse complement for guides occurring on rev strand
        guides_list = reverse_guides(guides_list, verbosity)
        # set variants for current guide
        guides_list = annotate_variants(guides_list, verbosity, debug)
        # annotate each guide with azimuth scores
        guides_list = azimuth_score(guides_list, verbosity, debug)
        # annotate each guide with rs3 scores
        guides_list = rs3_score(guides_list, verbosity, debug)
        # annotate each guide with CFDon scores
        guides_list = cfdon_score(guides_list, verbosity, debug)
        # annotate each guide functionally
        if functional_annotation:
            guides_list = funcann_guides(guides_list, region.contig, functional_annotation, "func", verbosity, debug)
        # annotate each guide with gene data
        if gene_annotation:
            guides_list = funcann_guides(guides_list, region.contig, gene_annotation, "gene", verbosity, debug)
        #
        #
        # if estimate_offtargets:  # estimate off-targets for each guide
        #     print("estimating off-targets")
        #     search_offtargets(guides_list, pam, genome, debug)
        #     exit()
        guides[region] = guides_list  # store annotated guides
    print_verbosity(
        f"Annotation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides
