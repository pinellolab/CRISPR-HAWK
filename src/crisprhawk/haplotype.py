""" """

from .crisprhawk_error import CrisprHawkHaplotypeError, CrisprHawkIupacTableError
from .exception_handlers import exception_handler
from .region import Region
from .sequence import Sequence
from .coordinate import Coordinate
from .variant import VariantRecord, VTYPES
from .utils import match_iupac, IUPAC_ENCODER

from typing import List, Dict, Tuple

import os


class Haplotype(Region):
    def __init__(
        self,
        sequence: Sequence,
        coord: Coordinate,
        phased: bool,
        chromcopy: int,
        debug: bool,
    ) -> None:
        self._debug = debug  # store debug flag
        super().__init__(sequence, coord)  # genomic sequence and region coordinates
        self._size = len(sequence)  # haplotype size
        self._variant_alleles = {}
        self._variants = "NA"  # haplotype variants
        self._samples = "REF"  # haplotype samples
        self._phased = phased  # haplotype phasing
        self._chromcopy = chromcopy  # chromosome copy
        self._initialize_posmap([], self._coordinates.start)  # initialize position map

    def __str__(self) -> str:
        return f"{self._samples}: {self._sequence.sequence}"

    def _initialize_posmap(self, chains: List[int], start: int) -> None:
        self._posmap = {
            p: start + i for i, p in enumerate(range(self._size + sum(chains)))
        }
        self._posmap_reverse = {pos: posrel for posrel, pos in self._posmap.items()}

    def _update_sequence(self, start: int, stop: int, alt: str) -> None:
        self._sequence._sequence_raw = (
            self._sequence._sequence_raw[:start]
            + list(alt.lower())
            + self._sequence._sequence_raw[stop:]
        )

    def substring(self, start: int, stop: int) -> str:
        return "".join(self._sequence._sequence_raw[start:stop])

    def _update_posmap(self, posrel: int, chain: int) -> None:
        if chain < 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] -= chain  # self._posmap[pos] + abs(chain)
        if chain > 0:
            for pos in range(posrel + 1, max(self._posmap.keys()) + 1):
                self._posmap[pos] = (
                    self._posmap[posrel]
                    if pos < posrel + chain + 1
                    else self._posmap[pos] - chain
                )
        self._posmap_reverse = {pos: posrel for posrel, pos in self._posmap.items()}

    def _update_variant_alleles(self, pos: int, offset: int) -> None:
        self._variant_alleles = {
            (p + offset if p > pos else p): alleles for p, alleles in self._variant_alleles.items()
        }    

    def _insert_variant_phased(
        self, position: int, ref: str, alt: str, chain: int, offset: int
    ) -> None:
        posrel = self._posmap_reverse[position]
        posrel_stop = posrel + abs(chain) + 1 if chain < 0 else posrel + 1
        if posrel_stop > self._size:
            posrel_stop = (self._size + offset) - 1
        refnt = self.substring(posrel, posrel_stop)
        if refnt != ref and refnt.isupper():
            raise ValueError(
                f"Mismatching reference alleles in VCF and reference sequence at position {position} ({refnt} - {ref})"
            )
        self._update_sequence(posrel, posrel_stop, alt)
        self._update_posmap(posrel, chain)

    def add_variants_phased(self, variants: List[VariantRecord], sample: str) -> None:
        if not self._phased:
            exception_handler(ValueError, "Unphased haplotype, unable to add phased variants", os.EX_DATAERR, True) # always trace this error
        variants = _sort_variants(variants)
        chains = _compute_chains(variants) # retrieve original sequence positions
        # position map to handle indels' presence
        self._initialize_posmap(chains, self._coordinates.start)
        for i, variant in enumerate(variants): # add variants to haplotype
            self._insert_variant_phased(
                variant.position,
                variant.ref,
                variant.alt[0],
                chains[i],
                sum(chains[:i]),
            )
        # add chromcopy as suffix to haplotype's samples
        suffix = "1|0" if self._chromcopy == 0 else "0|1"
        self._sequence = Sequence(
            "".join(self._sequence._sequence_raw), self._debug, allow_lower_case=True
        ) # reconstruct haplotype sequence
        self._samples = f"{sample}:{suffix}" if self._phased else sample
        self._variants = ",".join([v.id[0] for v in variants])

    def _insert_variant_unphased(
        self, position: int, ref: str, alt: str, vtype: str, chain: int, offset: int
    ):
        try:
            posrel = self._posmap_reverse[position]
        except KeyError:  # position may be deleted by previous deletion
            return
        posrel_stop = posrel + abs(chain) + 1 if chain < 0 else posrel + 1
        if posrel_stop > self._size:
            posrel_stop = (self._size + offset) - 1
        refnt = self.substring(posrel, posrel_stop) # retrieve ref nt in haplotype
        if not match_iupac(ref, refnt): # check ref alleles coherence
            raise ValueError(
                f"Mismatching reference alleles in VCF and reference sequence at position {position} ({refnt} - {ref})"
            )
        self._variant_alleles[posrel] = (ref, alt) # to solve haplotype
        if vtype == VTYPES[0]:  # if snv encode as iupac
            alt = _encode_iupac(ref, alt, position, self._debug)
        if vtype == VTYPES[1]:  # if indel update variant alleles map positions
            self._update_variant_alleles(posrel, len(alt) - len(ref))
        # update haplotype sequence and positions map
        self._update_sequence(posrel, posrel_stop, alt)
        self._update_posmap(posrel, chain)
        
    def add_variants_unphased(self, variants: List[VariantRecord], sample: str) -> None:
        variants = _sort_variants(variants)
        chains = _compute_chains(variants)
        self._initialize_posmap(chains, self._coordinates.start)
        for i, variant in enumerate(variants):
            self._insert_variant_unphased(
                variant.position,
                variant.ref,
                variant.alt[0],
                variant.vtype[0],
                chains[i],
                sum(chains[:i]),
            )
        self._sequence = Sequence(
            "".join(self._sequence._sequence_raw), self._debug, allow_lower_case=True
        )
        self._samples = sample
        self._variants = ",".join([v.id[0] for v in variants])

    def homozygous_samples(self) -> None:
        # if samples are homozygous, change their phasing value (support diploid)
        if self._samples == "REF":
            exception_handler(CrisprHawkHaplotypeError, "REF haplotype cannot be homozygous", os.EX_DATAERR, self._debug)  # type: ignore
        self._samples = ",".join(
            [
                f"{s[0]}:1|1"
                for sample in self._samples.split(",")
                for s in [sample.split(":")]
            ]
        )

    def set_samples(self, samples: str) -> None:
        self._samples = samples  # set samples to haplotype

    def set_variants(self, variants: str) -> None:
        self._variants = variants  # set variants to haplotype

    def set_posmap(self, posmap: Dict[int, int], posmap_rev: Dict[int, int]) -> None:
        self._posmap = posmap  # set position map to haplotype
        self._posmap_reverse = posmap_rev

    def set_id(self, hapid: str) -> None:
        self._id = hapid  # set haplotype ID

    def set_variant_alleles(self, variant_alleles: Dict[int, Tuple[str, str]]) -> None:
        self._variant_alleles = variant_alleles  # set variant alleles

    @property
    def samples(self) -> str:
        return self._samples

    @property
    def variants(self) -> str:
        return self._variants    

    @property
    def phased(self) -> bool:
        return self._phased

    @property
    def posmap(self) -> Dict[int, int]:
        return self._posmap

    @property
    def posmap_rev(self) -> Dict[int, int]:
        return self._posmap_reverse

    @property
    def id(self) -> str:
        return self._id
    
    @property
    def variant_alleles(self) -> Dict[int, Tuple[str, str]]:
        return self._variant_alleles


def _sort_variants(variants: List[VariantRecord]) -> List[VariantRecord]:
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


def _encode_iupac(ref: str, alt: str, position: int, debug: bool) -> str:
    try:
        return IUPAC_ENCODER["".join({ref, alt})]
    except KeyError as e:
        exception_handler(CrisprHawkIupacTableError, f"An error occurred while encoding {ref}>{alt} at position {position} as IUPAC character", os.EX_DATAERR, debug, e)  # type: ignore


class HaplotypeIndel(Haplotype):
    def __init__(self, sequence: Sequence, coord: Coordinate, phased: bool, chromcopy: int, debug: bool) -> None:
        super().__init__(sequence, coord, phased, chromcopy, debug)
        self._offset = 0
        self._indel_position = -1

    def set_offset(self, offset: int) -> None:
        self._offset = offset

    def set_indel_position(self, position: int) -> None:
        self._indel_position = position

    @property
    def offset(self) -> int:
        return self._offset
    
    @property
    def indel_position(self) -> int:
        return self._indel_position
