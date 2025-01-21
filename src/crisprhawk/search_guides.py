"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkBitsetError
from sequences import PAM
from bedfile import Region, IndelRegion
from bitset import Bitset

from typing import Tuple, List, Union

import os


def match(bitset1: Bitset, bitset2: Bitset) -> bool:
    # bitwise matching operation for input bitsets
    # assumes the two bitsets have the same length
    return all((ntbit & bitset2[i]).to_bool() for i, ntbit in enumerate(bitset1))


def scan_pam(
    pos: int,
    matches: List[int],
    bitset1: Bitset,
    bitset2: Bitset,
    strand: int,
    debug: bool,
) -> List[int]:
    try:  # scan pam occurrence at query position
        if match(bitset1, bitset2):
            matches.append(pos)  # pam occurrence found on current strand
    except ValueError as e:
        strand = "forward" if strand == 0 else "reverse"
        exception_handler(
            CrisprHawkBitsetError,
            f"PAM bitwise matching failed on {strand} strand",
            os.EX_DATAERR,
            debug,
            e,
        )
    return matches


def pam_search(
    pam: PAM, region: Union[Region, IndelRegion], guidelen: int, debug: bool
) -> Tuple[List[int], List[int]]:
    # assert that the input region has already bin encoded in bits
    assert hasattr(region, "_sequence_bits")
    pam.encode()  # encode pam linear time search
    # matching limit on the right - ignore pad nts
    stop_position = len(region) - len(pam) + 1
    # lists storing hits for input pam on forward and reverse strands
    matches_fwd, matches_rev = [], []
    # establish whether we scan a regular or an indel region
    indel = isinstance(region, IndelRegion)
    guidepamlen = guidelen + len(pam)  # used for indel region scan boundaries
    # start pam search from real region start - ignore pad nts
    # instead, if indel region start from the beginning
    scan_start = 0 if indel else guidelen
    for pos in range(scan_start, stop_position):
        for strand, matches, pam_bit in [
            (0, matches_fwd, pam.bits[0]),
            (1, matches_rev, pam.bits[1]),
        ]:
            if not indel or (
                indel
                and (
                    (strand == 0 and pos >= guidepamlen - 1)
                    or (strand == 1 and pos < guidepamlen)
                )
            ):
                scan_pam(
                    pos,
                    matches,
                    region.bits[pos : (pos + len(pam))],
                    pam_bit,
                    strand,
                    debug,
                )
    return matches_fwd, matches_rev


def extract_guide(
    region: Region, pos: int, guidelen: int, pamlen: int, right: bool
) -> str:
    if right:
        return region[pos : pos + guidelen + pamlen]
    return region[pos - guidelen : pos + pamlen]


def valid_position(pos: int, guidelen: int, regionlen: int, right: bool) -> bool:
    if right:
        return guidelen <= regionlen - pos
    return guidelen <= pos


def retrieve_guides(
    matches_fwd: List[int],
    matches_rev: List[int],
    pamlen: int,
    region: Region,
    guidelen: int,
    right: bool,
) -> Tuple[List[List[str]], List[List[str]]]:
    # helper function to recover guides sequence
    def get_guides(matches: List[int], guideright: bool) -> List[List[str]]:
        return [
            extract_guide(region, pos, guidelen, pamlen, guideright)
            for pos in matches
            if valid_position(pos, guidelen, len(region), guideright)
        ]

    guides_fwd = get_guides(matches_fwd, right)  # guides on positive strand
    guides_rev = get_guides(matches_rev, not right)  # guides on negative strand
    return guides_fwd, guides_rev


def search(
    pam: PAM, region: Region, guidelen: int, right: bool, debug: bool
) -> Tuple[Tuple[List[int], List[int]], Tuple[List[List[str]], List[List[str]]]]:
    # search pam occurrences on forward and reverse strand of the input sequence
    matches_fwd, matches_rev = pam_search(pam, region, guidelen, debug)
    # recover guide candidates found on forward and reverse strands
    guides_fwd, guides_rev = retrieve_guides(
        matches_fwd, matches_rev, len(pam), region, guidelen, right
    )
    return (matches_fwd, matches_rev), (guides_fwd, guides_rev)
