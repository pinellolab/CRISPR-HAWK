"""This module provides functions for searching and extracting CRISPR guide
sequences from genomic regions and haplotypes.

It includes utilities for scanning for PAM occurrences, extracting and validating
guide sequences, handling IUPAC ambiguity codes, and removing redundant guides.
The module supports both reference and alternative haplotypes, phased and unphased
data, and is designed to facilitate flexible and robust guide discovery for genome
editing applications.
"""

from .exception_handlers import exception_handler
from .crisprhawk_error import (
    CrisprHawkBitsetError,
    CrisprHawkIupacTableError,
    CrisprHawkCfdScoreError,
)
from .utils import print_verbosity, flatten_list, VERBOSITYLVL, IUPACTABLE
from .guide import Guide, GUIDESEQPAD
from .bitset import Bitset
from .pam import PAM
from .region_constructor import PADDING
from .region import Region
from .haplotype import Haplotype, HaplotypeIndel

from collections import defaultdict

from typing import List, Tuple, Dict, Union, Optional, DefaultDict
from itertools import product
from time import time

import os


def match(
    bitset1: List[Bitset], bitset2: List[Bitset], position: int, debug: bool
) -> bool:
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
            CrisprHawkBitsetError,  # type: ignore
            f"PAM bitwise matching failed at position {position}",
            os.EX_DATAERR,
            debug,
            e,
        )


def compute_scan_start_stop(
    hap: Haplotype, region_start: int, region_stop: int, pamlen: int
) -> Tuple[int, int]:
    """Computes the scan start and stop positions for a haplotype within a genomic
    region.

    Determines the sequence coordinates to scan for PAMs, accounting for region
    padding and indels at region boundaries.

    Args:
        hap (Haplotype): The haplotype object containing position mappings.
        region_start (int): The start coordinate of the region.
        region_stop (int): The stop coordinate of the region.
        pamlen (int): The length of the PAM sequence.

    Returns:
        Tuple[int, int]: The start and stop positions for scanning in the haplotype
            sequence.
    """
    stop_p = min(region_stop - PADDING, hap.stop)  # stop position with padding
    # handle indels overlapping the end of region
    if stop_p == region_stop - PADDING and stop_p not in hap.posmap_rev:
        # find the highest available position in the mapping
        upperbound_p = max(hap.posmap_rev.keys())
        # search for next available position in the mapping
        search_range = upperbound_p - stop_p + 1
        for offset in range(search_range):
            p = stop_p + offset
            if p in hap.posmap_rev:
                stop_p = p
                break
    # convert genomic coordinates to sequence coordinates
    scan_stop = hap.posmap_rev[stop_p] - pamlen + 1
    start_p = max(region_start + PADDING, hap.start)
    scan_start = hap.posmap_rev[start_p]
    return scan_start, scan_stop


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
        scan_start, scan_stop = compute_scan_start_stop(
            hap, region.start, region.stop, len(pam)
        )
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


def _valid_guide(
    pamguide: str, pam: PAM, direction: int, right: bool, debug: bool
) -> bool:
    """Checks if a guide sequence is valid for a given PAM and direction.

    Encodes the guide as a PAM and verifies if it matches the PAM bits for the
    specified strand direction.

    Args:
        pamguide (str): The guide sequence to validate.
        pam (PAM): The PAM object to match against.
        direction (int): The strand direction (0 for positive, 1 for negative).
        right (bool): Whether the guide is on the right strand.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        bool: True if the guide is valid for the given PAM and direction, False
            otherwise.
    """
    p = PAM(pamguide, right, debug)
    p.encode(0)
    if direction == 0:  # positive
        return all((ntbit & pam.bits[i]).to_bool() for i, ntbit in enumerate(p.bits))
    return all((ntbit & pam.bitsrc[i]).to_bool() for i, ntbit in enumerate(p.bits))


def _decode_iupac(nt: str, pos: int, llimit: int, rlimit: int, h: Haplotype, debug: bool) -> str:
    """Decodes an IUPAC nucleotide code at a given position within specified limits.

    Returns the resolved nucleotide(s) for the given IUPAC code, handling ambiguity 
    and variant alleles, or the original nucleotide if out of bounds.

    Args:
        nt (str): The nucleotide or IUPAC code to decode.
        pos (int): The position of the nucleotide in the sequence.
        llimit (int): The left boundary for decoding.
        rlimit (int): The right boundary for decoding.
        h (Haplotype): The Haplotype object containing variant allele information.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        str: The decoded nucleotide(s) or the original nucleotide if out of bounds.

    Raises:
        CrisprHawkIupacTableError: If the IUPAC character is invalid.
    """
    if pos < llimit or pos > rlimit:
        return nt
    try:
        ntiupac = IUPACTABLE[nt.upper()]
    except KeyError as e:
        exception_handler(CrisprHawkIupacTableError, f"Invalid IUPAC character ({nt})", os.EX_DATAERR, debug, e)  # type: ignore
    if len(ntiupac) == 1:  # A, C, G, or T (reference / indel / homozygous)
        return ntiupac.lower() if nt.islower() else ntiupac
    return "".join(
        [n if n == h.variant_alleles[pos][0] else n.lower() for n in list(ntiupac)]
    )


def resolve_guide(
    guideseq: str,
    pam: PAM,
    direction: int,
    right: bool,
    pos: int,
    guidelen: int,
    h: Haplotype,
    debug: bool,
) -> List[str]:
    """Resolves all possible guide sequence alleles for a given IUPAC-encoded guide.

    Decodes IUPAC ambiguity codes in the guide sequence to generate all possible
    allele combinations, filtering for valid PAM matches.

    Args:
        guideseq (str): The IUPAC-encoded guide sequence.
        pam (PAM): The PAM object to match against.
        direction (int): The strand direction (0 for positive, 1 for negative).
        right (bool): Whether the guide is on the right strand.
        pos (int): The position in the haplotype.
        guidelen (int): The length of the guide sequence.
        h (Haplotype): The Haplotype object containing variant allele information.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        List[str]: A list of all valid guide sequence alleles matching the PAM.
    """
    # retrieve guide sequence relative start position
    p = pos - GUIDESEQPAD if right else pos - guidelen - GUIDESEQPAD
    llimit = pos - 10 if right else pos - guidelen - 10
    rlimit = pos + len(pam) + guidelen + 10 if right else pos + len(pam) + 10 
    guide_alts = [
        "".join(g)
        for g in product(
            *[list(_decode_iupac(nt, p + i, llimit, rlimit, h, debug)) for i, nt in enumerate(guideseq)]
        )
    ]  # decode iupac string
    idx = GUIDESEQPAD if right else (len(guideseq) - GUIDESEQPAD - len(pam))
    return [
        guide
        for guide in guide_alts
        if _valid_guide(guide[idx : idx + len(pam)], pam, direction, right, debug)
    ]


def adjust_guide_position(
    posmap: Dict[int, int], posrel: int, guidelen: int, pamlen: int, right: bool
) -> Tuple[int, int]:
    """Adjusts and returns the genomic start and stop positions for a guide sequence.

    Calculates the start and stop positions in the genome based on the relative
    position, guide length, PAM length, and strand orientation.

    Args:
        posmap (Dict[int, int]): Mapping from sequence-relative to genomic positions.
        posrel (int): The relative position in the sequence.
        guidelen (int): The length of the guide sequence.
        pamlen (int): The length of the PAM sequence.
        right (bool): Whether the guide is on the right (forward) strand.

    Returns:
        Tuple[int, int]: The genomic start and stop positions for the guide.
    """
    start = posmap[posrel] if right else posmap[posrel - guidelen]
    stop = posmap[posrel + guidelen + pamlen] if right else posmap[posrel + pamlen]
    return start, stop


def group_guides_position(
    guides: List[Guide], debug: bool
) -> DefaultDict[str, Dict[int, Union[Optional[Guide], List[Guide]]]]:
    """Groups guides by their genomic start position and strand.

    Organizes guides into a dictionary keyed by position and strand, separating
    reference and alternative guides for downstream analysis.

    Args:
        guides (List[Guide]): A list of Guide objects to group.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        DefaultDict[str, Dict[int, Union[Optional[Guide], List[Guide]]]]:
            A dictionary mapping position/strand keys to reference and
            alternative guides.
    """
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


def remove_redundant_guides(guides: List[Guide], debug: bool) -> List[Guide]:
    """Removes redundant guide sequences from a list of guides based on sequence
    and sample origin.

    Groups guides by position and strand, then filters out guides that are redundant
    with the reference guide at each position.

    Args:
        guides (List[Guide]): A list of Guide objects to filter.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        List[Guide]: A list of non-redundant Guide objects.
    """
    filtered_guides = []  # list containing non-redundant guides
    grouped_guides = group_guides_position(guides, debug)  # group guides by position
    for pos, guide_group in grouped_guides.items():
        refguide, altguides = guide_group[0], guide_group[1]
        if refguide is None:  # only alt guides, redundancy not possible
            filtered_guides.extend(altguides)  # type: ignore
            continue
        altguides_list: List[Guide] = altguides  # type: ignore
        for guide in altguides_list:  # remove redundant guides
            refguide_seq = refguide.sequence[GUIDESEQPAD:-GUIDESEQPAD].upper()  # type: ignore
            guide_seq = guide.sequence[GUIDESEQPAD:-GUIDESEQPAD].upper()
            if (guide.samples != "REF" and guide_seq != refguide_seq) or (
                guide.samples == "REF" and guide_seq == refguide_seq
            ):
                filtered_guides.append(guide)
    return filtered_guides


def is_pamhit_valid(
    pamhit_pos: int, haplen: int, guidelen: int, pamlen: int, right: bool
) -> bool:
    """Checks if a PAM hit position is valid for guide extraction within a haplotype.

    Determines whether the guide sequence can be safely extracted from the haplotype
    without exceeding sequence boundaries.

    Args:
        pamhit_pos (int): The position of the PAM hit in the haplotype.
        haplen (int): The length of the haplotype sequence.
        guidelen (int): The length of the guide sequence.
        pamlen (int): The length of the PAM sequence.
        right (bool): Whether the guide is on the right (forward) strand.

    Returns:
        bool: True if the PAM hit position is valid for guide extraction, False otherwise.
    """
    if right:
        return pamhit_pos + guidelen + pamlen + GUIDESEQPAD < haplen
    return pamhit_pos - guidelen - GUIDESEQPAD >= 0


def is_pamhit_in_range(
    poshit: int, guidelen: int, pamlen: int, haplen: int, right: bool
) -> bool:
    """Checks if a PAM hit position is within the valid range for guide extraction.

    Computes the left and right boundaries for the guide sequence and ensures they
    are within the haplotype length.

    Args:
        poshit (int): The position of the PAM hit in the haplotype.
        guidelen (int): The length of the guide sequence.
        pamlen (int): The length of the PAM sequence.
        haplen (int): The length of the haplotype sequence.
        right (bool): Whether the guide is on the right (forward) strand.

    Returns:
        bool: True if the guide sequence is within the valid range, False otherwise.
    """
    # compute left and right boundaries for current guide
    lbound = poshit - GUIDESEQPAD if right else poshit - guidelen - GUIDESEQPAD
    rbound = (
        poshit + guidelen + pamlen + GUIDESEQPAD
        if right
        else poshit + pamlen + GUIDESEQPAD
    )
    return lbound >= 0 and rbound <= haplen  # check whether guide sequence is in range


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
    """Retrieves guide sequences from a haplotype at specified PAM hit positions.

    Extracts and returns a list of Guide objects for each valid PAM hit, handling
    IUPAC ambiguity and filtering for reference and alternative guides as needed.

    Args:
        pam_hits (List[int]): List of PAM hit positions in the haplotype.
        haplotype (Haplotype): The Haplotype object to extract guides from.
        guidelen (int): The length of the guide sequence.
        pam (PAM): The PAM object to match against.
        direction (int): The strand direction (0 for positive, 1 for negative).
        right (bool): Whether the guide is on the right (forward) strand.
        variants_present (bool): Whether variants are present in the haplotype.
        phased (bool): Whether the haplotype is phased.
        verbosity (int): Verbosity level for logging.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        List[Guide]: A list of Guide objects extracted from the haplotype.
    """
    guides = []  # list of haplotype-specific guides
    s = "-" if direction == 1 else "+"
    print_verbosity(
        f"Retrieving guides from haplotype {haplotype.samples} (strand: {s})",
        verbosity,
        VERBOSITYLVL[3],
    )
    start = time()
    for pos in pam_hits:
        if not is_pamhit_in_range(pos, guidelen, len(pam), len(haplotype), right):
            continue  # if guide sequence not in range skip the hit
        # retrieve guide sequence from haplotype
        guideseq = extract_guide_sequence(haplotype, pos, len(pam), guidelen, right)
        if (
            haplotype.samples != "REF" and guideseq[GUIDESEQPAD:-GUIDESEQPAD].isupper()
        ):  # reference guide
            continue
        guideseqs = [guideseq]
        if variants_present and not phased:  # handle guides with iupac
            if is_pamhit_valid(pos, len(haplotype), guidelen, len(pam), right):
                guideseqs = resolve_guide(
                    guideseq, pam, direction, right, pos, guidelen, haplotype, debug
                )
            else:
                guideseqs = []
        for guideseq in guideseqs:
            # compute guide's start and stop positions
            guide_start, guide_stop = adjust_guide_position(
                haplotype.posmap, pos, guidelen, len(pam), right
            )
            guide = Guide(
                guide_start,
                guide_stop,
                guideseq,
                guidelen,
                len(pam),
                direction,
                haplotype.samples,
                haplotype.variants,
                haplotype.afs,
                debug,
                right,
                haplotype.id,
            )
            guides.append(guide)  # report guide
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
    """Searches for guide candidates in a genomic region across multiple haplotypes.

    Identifies all valid guide sequences matching the PAM in the specified region
    and haplotypes, returning a list of non-redundant Guide objects.

    Args:
        pam (PAM): The PAM object to search for.
        region (Region): The genomic region to search within.
        haplotypes (List[Haplotype]): List of Haplotype objects to scan.
        haplotypes_bits (List[List[Bitset]]): Bit-encoded haplotype sequences.
        guidelen (int): The length of the guide sequence.
        right (bool): Whether to search on the right (forward) strand.
        variants_present (bool): Whether variants are present in the haplotypes.
        phased (bool): Whether the haplotypes are phased.
        verbosity (int): Verbosity level for logging.
        debug (bool): Whether to enable debug mode for error handling.

    Returns:
        List[Guide]: A list of non-redundant Guide objects found in the region.
    """
    print_verbosity(
        f"Searching guide candidates in {region.coordinates}",
        verbosity,
        VERBOSITYLVL[3],
    )
    # search pam occurrences on forward and reverse strand of the input sequence
    pam_hits = pam_search(pam, region, haplotypes, haplotypes_bits, verbosity, debug)
    # retrieve guide candidates sequence for each haplotype
    guides = flatten_list(
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
    return remove_redundant_guides(guides, debug)  # remove redundant guides
