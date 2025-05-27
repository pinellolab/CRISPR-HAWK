""" """

from .exception_handlers import exception_handler
from .crisprhawk_error import CrisprHawkBitsetError, CrisprHawkIupacTableError
from .utils import print_verbosity, flatten_list, VERBOSITYLVL, IUPACTABLE
from .guide import Guide, GUIDESEQPAD
from .bitset import Bitset
from .pam import PAM
from .region_constructor import PADDING
from .region import Region
from .haplotype import Haplotype

from typing import List, Tuple, Dict
from itertools import product
from time import time

import os


def match(bitset1: List[Bitset], bitset2: List[Bitset], position: int, debug: bool) -> bool:
    """Performs a bitwise match between two Bitset objects at a given position.

    This function checks if all corresponding bits match between two bitsets, 
    raising an exception if the operation fails.

    Args:
        bitset1: The first Bitset to compare.
        bitset2: The second Bitset to compare.
        position: The position in the sequence for error reporting.
        debug: Whether to enable debug mode for error handling.

    Returns:
        bool: True if all bits match, False otherwise.

    Raises:
        CrisprHawkBitsetError: If the bitwise matching operation fails.
    """
    # bitwise matching operation for input bitsets
    # assumes the two bitsets have the same length
    try:
        return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))
    except ValueError as e:
        exception_handler(
            CrisprHawkBitsetError, # type: ignore
            f"PAM bitwise matching failed at position {position}",
            os.EX_DATAERR,
            debug,
            e,
        )


def scan_haplotype(
    pam: PAM, haplotype: List[Bitset], start: int, stop: int, debug: bool
) -> Tuple[List[int], List[int]]:
    """Scans a haplotype for PAM matches on both forward and reverse strands.

    This function returns the positions of all matches for the PAM and its 
    reverse complement within the specified range.

    Args:
        pam: The PAM object containing bit representations for matching.
        haplotype: The haplotype sequence as a list of Bitset objects.
        start: The start position for scanning.
        stop: The stop position for scanning.
        debug: Whether to enable debug mode for error handling.

    Returns:
        Tuple[List[int], List[int]]: Lists of match positions for forward and 
            reverse strands.
    """
    # lists storing hits for input pam on forward and reverse strands
    matches_fwd, matches_rev = [], []
    for pos in range(start, stop):  # start scanning haplotype for guides
        if match(pam.bits, haplotype[pos : (pos + len(pam))], pos, debug):
            matches_fwd.append(pos)  # hit on forward strand
        if match(pam.bitsrc, haplotype[pos : (pos + len(pam))], pos, debug):
            matches_rev.append(pos)  # hit on negative strand
    return matches_fwd, matches_rev


def pam_search(
    pam: PAM,
    region: Region,
    haplotypes: List[Haplotype],
    haplotypes_bits: List[List[Bitset]],
    verbosity: int,
    debug: bool,
) -> List[Tuple[List[int], List[int]]]:
    """Searches for PAM occurrences in each haplotype within a genomic region.

    This function scans all haplotypes for matches to the PAM and its reverse 
    complement, returning the positions of all matches.

    Args:
        pam: The PAM object containing bit representations for matching.
        region: The genomic region to search within.
        haplotypes: A list of Haplotype objects to scan.
        haplotypes_bits: A list of bit-encoded haplotype sequences.
        verbosity: The verbosity level for logging.
        debug: Whether to enable debug mode for error handling.

    Returns:
        List[Tuple[List[int], List[int]]]: A list of tuples containing forward and 
            reverse strand match positions for each haplotype.
    """
    pam_hits = []  # list of pam hits
    for i, hap in enumerate(haplotypes):
        print_verbosity(
            f"Searching PAM occurrences in haplotype {hap.samples}",
            verbosity,
            VERBOSITYLVL[3],
        )
        # define scan stop position for each haplotype
        scan_stop = hap.posmap_rev[region.stop - PADDING] - len(pam) + 1
        scan_start = hap.posmap_rev[region.start + PADDING]  # start position for guides search
        # scan haplotype for pam occurrences
        matches_fwd, matches_rev = scan_haplotype(
            pam, haplotypes_bits[i], scan_start, scan_stop, debug
        )
        print_verbosity(
            f"Found {len(matches_fwd) + len(matches_rev)} PAM occurrences ({len(matches_fwd)} on 5'-3'; {len(matches_rev)} on 3'-5')",
            verbosity,
            VERBOSITYLVL[3],
        )
        pam_hits.append((matches_fwd, matches_rev))  # store pam hits for each haplotype
    return pam_hits


def extract_guide_sequence(
    haplotype: Haplotype, position: int, pamlen: int, guidelen: int, right: bool
) -> str:
    """Extracts the guide RNA sequence from a haplotype at a specified position.

    This function returns the guide sequence including the PAM and optional 
    padding, depending on the strand direction.

    Args:
        haplotype: The Haplotype object to extract the guide from.
        position: The position in the haplotype to start extraction.
        pamlen: The length of the PAM sequence.
        guidelen: The length of the guide sequence.
        right: Whether to extract from the right (forward) strand.

    Returns:
        str: The extracted guide sequence.
    """
    if right:
        return "".join(
            haplotype[
                position - GUIDESEQPAD : position + guidelen + pamlen + GUIDESEQPAD
            ]
        )
    return "".join(
        haplotype[position - guidelen - GUIDESEQPAD : position + pamlen + GUIDESEQPAD]
    )

def _valid_guide(pamguide: str, pam: PAM, direction: int, debug: bool) -> bool:
    p = PAM(pamguide, debug)
    p.encode(0)
    if direction == 0:  # positive
        return all((ntbit & pam.bits[i]).to_bool() for i, ntbit in enumerate(p.bits))
    return all((ntbit & pam.bitsrc[i]).to_bool() for i, ntbit in enumerate(p.bits))

def _decode_iupac(nt: str, debug: bool) -> str:
    try:
        ntiupac = IUPACTABLE[nt.upper()]
    except KeyError as e:
        exception_handler(CrisprHawkIupacTableError, f"Invalid IUPAC character ({nt})", os.EX_DATAERR, debug, e) # type: ignore
    return ntiupac.lower() if nt.islower() else ntiupac

def resolve_guide(guideseq: str, pam: PAM, direction: int, right: bool, debug: bool) -> List[str]:
    guide_alts = ["".join(g) for g in product(*[list(_decode_iupac(nt, debug)) for nt in guideseq])]
    idx = GUIDESEQPAD if right else (len(guideseq) - GUIDESEQPAD - len(pam))
    return [guide for guide in guide_alts if _valid_guide(guide[idx:idx + len(pam)], pam, direction, debug)]


def adjust_guide_position(posmap: Dict[int, int], posrel: int, guidelen: int, pamlen: int, right: bool) -> Tuple[int, int]:
    start = posmap[posrel] if right else posmap[posrel - guidelen]
    stop = posmap[posrel + guidelen + pamlen] + 1 if right else posmap[posrel + pamlen]
    return start, stop

def retrieve_guides(
    pam_hits: List[int],
    haplotype: Haplotype,
    guidelen: int,
    pam: PAM,
    direction: int,
    right: bool,
    variants_present: bool,
    phased: bool,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    guides = []  # list of haplotype-specific guides
    print_verbosity(
        f"Retrieving guides from haplotype {haplotype.samples}",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()
    for pos in pam_hits:
        # retrieve guide sequence from haplotype
        guideseq = extract_guide_sequence(haplotype, pos, len(pam), guidelen, right)
        if (
            haplotype.samples != "REF" and guideseq[GUIDESEQPAD:-GUIDESEQPAD].isupper()
        ):  # reference guide
            continue
        guideseqs = [guideseq]
        if variants_present and not phased:  # handle guides with iupac 
            guideseqs = resolve_guide(guideseq, pam, direction, right, debug)
        for guideseq in guideseqs:
            # compute guide's start and stop positions
            guide_start, guide_stop = adjust_guide_position(haplotype.posmap, pos, guidelen, len(pam), right)
            guides.append(
                Guide(
                    guide_start,
                    guide_stop,
                    guideseq,
                    guidelen,
                    len(pam),
                    direction,
                    haplotype.samples,
                    haplotype.variants,
                    debug,
                    right,
                    haplotype.id,
                )
            )
    print_verbosity(
        f"Guides retrieved in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def search(
    pam: PAM,
    region: Region,
    haplotypes: List[Haplotype],
    haplotypes_bits: List[List[Bitset]],
    guidelen: int,
    right: bool,
    variants_present: bool,
    phased: bool,
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    print_verbosity(
        f"Searching guide candidates in {region.coordinates}",
        verbosity,
        VERBOSITYLVL[3],
    )
    # search pam occurrences on forward and reverse strand of the input sequence
    pam_hits = pam_search(pam, region, haplotypes, haplotypes_bits, verbosity, debug)
    # retrieve guide candidates sequence for each haplotype
    return flatten_list(
        [
            retrieve_guides(
                pam_hits[i][strand],
                haplotype,
                guidelen,
                pam,
                strand,
                (not right if strand == 1 else right),
                variants_present,
                phased,
                verbosity,
                debug,
            )
            for i, haplotype in enumerate(haplotypes)
            for strand in [0, 1]
        ]
    )
