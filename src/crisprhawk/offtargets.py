""" """

from .crisprhawk_error import CrisprHawkOffTargetsError
from .exception_handlers import exception_handler
from .utils import suppress_stdout, suppress_stderr
from .guide import Guide
from .pam import pam

from pybedtools import BedTool

import tempfile
import subprocess
import os

MAXMM = 4  # maximum number of mismatches allowed when searching offtargets
MAXENTRIES = 2000000  # maximum entries (bwa param)
MAXOCC = 60000  # maximum number of occurrences (bwa param)

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

def _read_sam_line(line: str, debug: bool) -> str:
    fields = line.strip().split()  # retrieve individual fields
    if fields[5] == "*":
        return ""  # skip unmapped (empty cigar)
    # retrieve mm, alignment score, alt alignments
    mm, alnscore, hasalt = _parse_tags(fields[11:])
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


def _sam2bed(guide_sam: str, guide_pos: int, pamlen: int, guidelen: int, right: bool, debug: bool) -> str:
    guide_bed = f"{guide_sam}.bed"
    try:  # retrieve bed file for found off-targets from bwa sam
        with open(guide_sam, mode="r") as infile, open(guide_bed, mode="w") as outfile:
            for line in infile:
                if line.startswith("@"):  # skip comments (sam)
                    continue
                if not (bedline := _read_sam_line(line, debug)):
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
    guide_sam = _samse(guide_fa, guide_sa, genome, debug=False)
    # generate bed file from sam
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


def retrieve_offtargets(guide_bed: str, genome: str, pam: PAM, pamlen: int, right: bool) -> List[Tuple[str, str, str, str, str]]:
    offtargets = offtarget_sequence(guide_bed, genome)  # retrieve offtargets sequence
    return filter_offtargets(offtargets, pam, pamlen, right)

def search_offtargets(guides: List[Guide], pam: PAM, genome: str, debug: bool):
    workingdir = tempfile.mkdtemp(prefix="crisprhawk_")  # create temp dir
    print(workingdir)
    try:
        for guide in guides:
            guide_fa = create_guide_fa(guide, workingdir, debug)  # guide fasta
            guide_sa = search(guide_fa, genome, guide.guidelen, debug)  # search offtargets
            guide_bed = sam2bed(guide_fa, guide_sa, genome, guide.start, guide.pamlen, guide.guidelen, guide.right, debug)  # convert sam to bed
            x = retrieve_offtargets(guide_bed, genome, pam, guide.pamlen, guide.right)
            print(len(x))
            print("\n".join(["\t".join(e) for e in x]))
            break

    except Exception:
        subprocess.call(f"rm -r {workingdir}", shell=True)  # remove temp dir
        raise Exception("Off-target search failed")
