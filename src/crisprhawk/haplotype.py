""" """

from region import Region
from sequence import Sequence
from coordinate import Coordinate
from variant import VariantRecord, VTYPES

from typing import List, Set


class Haplotype(Region):
    def __init__(self, sequence: Sequence, coord: Coordinate, phased: bool, chromcopy: int) -> None:
        super().__init__(sequence, coord)  # genomic sequence and region coordinates
        self._size = len(sequence)  # haplotype size
        self._variants = "NA"  # haplotype variants
        self._samples = "REF"  # haplotype samples
        self._phased = phased  #haplotype phasing
        self._chromcopy = chromcopy  # chromosome copy

    
    def _initialize_posmap(self, chains: List[int], start: int) -> None:
        self._posmap = {p: start + i for i, p in enumerate(range(sum(chains)))}


    def _update_sequence(self, start: int, stop: int, alt: str) -> None:
        self._sequence._sequence_raw = self._sequence._sequence_raw[:start] + list(alt.lower()) + self._sequence._sequence_raw[stop:]

    def substring(self, start: int, stop: int) -> str:
        return "".join(self._sequence._sequence_raw[start:stop]) 

    def _insert_variant(self, position: int, ref: str, alt: str, chain: int, offset: int):
        posrel = _adjust_position(position - 1, self._coordinates.start, offset)
        posrel_stop = posrel + abs(chain) + 1 if chain < 0 else posrel + 1
        if posrel_stop > self._size:
            posrel_stop = (self._size + offset) - 1
        refnt = self.substring(posrel, posrel_stop)
        if refnt != ref:
            raise ValueError(f"Mismatching reference alleles in VCF and reference sequence at position {position} ({refnt} - {ref})")
        self._update_sequence(posrel, posrel_stop, alt)
        self._update_posmap(posrel, position, chain)

    
    def _update_posmap(self, posrel: int, position: int, chain: int):
        if chain < 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] += chain
        if chain > 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] = self._posmap[posrel] if pos < posrel + chain + 1 else self._posmap[pos] - chain + 1


    
    def add_variants(self, variants: List[VariantRecord], sample: str):
        variants = _sort_variants(variants)
        chains = _compute_chains(variants)
        self._initialize_posmap(chains, self._coordinates.start)
        for i, variant in enumerate(variants):
            self._insert_variant(variant.position, variant.ref, variant.alt[0], chains[i], sum(chains[:i]))
        suffix = "" 
        if self._phased:
            suffix = "1|0" if self._chromcopy == 0 else "0|1"
        self._samples = f"{sample}:{suffix}" if self._phased else sample
        self._variants = ",".join([v.id[0] for v in variants])
        self.



def _sort_variants(variants: Set[VariantRecord]) -> List[VariantRecord]:
    # sort variants set to have snps before indels
    snps, indels = [], []
    for variant in variants:
        if variant.vtype[0] == VTYPES[0]:  # snp
            snps.append(variant)
        else:  # indel
            indels.append(variant)
    return sorted(snps) + sorted(indels)


def _compute_chains(variants: List[VariantRecord]) -> List[int]:
    return [len(v.alt[0]) - len(v.ref) for v in variants]

def _adjust_position(position: int, pivot_pos: int, offset: int) -> int:
    return (position - pivot_pos) + offset

def _compute_indel_length(ref: str, alt: str) -> int:
    """Computes the length of an indel.

    Args:
        ref: The reference sequence.
        alt: The alternate sequence.
    Returns:
        The length of the indel.
    """
    return abs(len(ref) - len(alt))
