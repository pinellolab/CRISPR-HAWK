"""
"""

from crisprhawk_argparse import CrisprHawkInputArgs
from regions import construct_regions
from haplotypes import reconstruct_haplotypes, reconstruct_haplotypes_ref
from utils import print_verbosity, VERBOSITYLVL
from search_guides import search
from annotate import annotate_guides
from encoder import encode
from bitset import Bitset
from guide import Guide
from pam import PAM

from hapsolver import Region, Haplotype
from typing import Dict, List
from time import time


def encode_pam(pam: str, verbosity: int, debug: bool) -> PAM:
    # construct pam object
    print_verbosity(f"Creating PAM object for PAM {pam}", verbosity, VERBOSITYLVL[1])
    start = time()  # encoding start time
    pam = PAM(pam, debug)
    pam.encode(verbosity)  # encode pam sequence
    print_verbosity(f"PAM object for PAM {pam} created in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2])
    return pam

def encode_haplotypes(haplotypes: Dict[Region, List[Haplotype]], verbosity: int, debug: bool) -> Dict[Region, List[List[Bitset]]]:
    # encode haplotypes in bit for efficient guide search
    print_verbosity("Encoding haplotypes in bits", verbosity, VERBOSITYLVL[1])
    start = time()  # encoding start time
    haplotypes_bits = {
        region: [encode(hap.sequence, verbosity, debug) for hap in haps]
        for region, haps in haplotypes.items()
    }  # encode input haplotypes as sequences of bits
    print_verbosity(f"Haplotype encoding completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2])
    return haplotypes_bits

def guides_search(pam: PAM, haplotypes: Dict[Region, List[Haplotype]], haplotypes_bits: Dict[Region, List[List[Bitset]]], guidelen: int, right: bool, verbosity: int, debug: bool) -> Dict[Region, List[Guide]]:
    # search guide candidates on encoded haplotypes
    print_verbosity("Searching guides on haplotypes", verbosity, VERBOSITYLVL[1])
    start = time()  # search start time
    guides = {
        region: search(pam, region, haplotype, haplotypes_bits[region], guidelen, right, verbosity, debug)
        for region, haplotype in haplotypes.items()
    }
    print_verbosity(f"Guides search completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2])
    return guides

def crisprhawk(args: CrisprHawkInputArgs) -> None:
    # extract genomic regions defined in input bed file
    regions = construct_regions(args.fasta, args.bedfile, args.fasta_idx, args.guidelen, args.verbosity, args.debug)
    if args.vcfs:  # establish whether variants have been given
        # reconstruct haplotypes in each region
        haplotypes = reconstruct_haplotypes(args.vcfs, regions, args.verbosity, args.debug)
    else:
        # reconstruct haplotypes with reference sequence only
        haplotypes = reconstruct_haplotypes_ref(regions, args.verbosity, args.debug)
    # encode pam and haplotype sequences in bit for efficient guides search
    pam = encode_pam(args.pam, args.verbosity, args.debug)
    haplotypes_bits = encode_haplotypes(haplotypes, args.verbosity, args.debug)
    # search guide candidates within input regions
    guides = guides_search(pam, haplotypes, haplotypes_bits, args.guidelen, args.right, args.verbosity, args.debug)
    # annotate guide candidates within each region
    guides = annotate_guides(guides, args.debug)
