"""
"""

from exception_handlers import exception_handler
from crisprhawk_error import CrisprHawkHaplotypeError
from variants import VariantRecord, VTYPES, INDELTYPES
from bedfile import IndelRegion, Region
from variant_map import VariantMap
from utils import STRAND
from guide import Guide

from typing import Tuple, List, Dict, Set, Union, Optional

import os

# iupac alphabet characters excluded canonical nucleotides [A, C, G, T, N]
IUPACSPECIALCHARS = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D"]


def find_variant_pos_indel(
    guide: str,
    position: int,
    variants_maps: Dict[int, VariantRecord],
    indelregion: IndelRegion,
    guidelen: int,
    pamlen: int,
) -> Tuple[str, Dict[int, VariantRecord]]:
    # store variants found within the guide; variants are denoted by iupac chars
    # different from the canonical nucleotides
    variants = {
        i: variants_maps[position + i]
        for i, _ in enumerate(guide)
        if position + i in variants_maps
    }
    guidepamlen = guidelen + pamlen  # total guide + pam length
    # compute guide reference sequence
    offset = (
        position if position <= len(indelregion) - guidepamlen else position - guidelen
    )
    indelposrel = indelregion.indelpos - offset
    indelposrel = 0 if indelposrel < 0 else indelposrel

    refseq = []
    i = 0
    while i < len(guide):
        if i in variants:
            if i == indelposrel:
                refseq.append(indelregion._ref)
                if indelregion.indel_type == INDELTYPES[0]:
                    i += len(indelregion._alt)
                    continue
            else:
                refseq.append(variants[i].ref)
        else:
            refseq.append(guide[i])
        i += 1
    refseq = "".join(refseq)
    return refseq, variants


def filter_guides_indel(
    sample_guides_chromcopy: Dict[str, Tuple[List[str], List[str], List[str]]]
) -> Dict[str, Tuple[List[str], List[str]]]:
    return {
        sample: (guide, list(dict.fromkeys(variants)))
        for sample, (guide, variants, vtypes) in sample_guides_chromcopy.items()
        if variants and VTYPES[1] in vtypes
    }


def _compute_indel_posrel(
    indelpos: int, region_len: int, guidelen: int, pamlen: int, guidepos: int
) -> int:
    guidepamlen = guidelen + pamlen  # total guide + pam length
    # adjust subtraction term to identify indel position, based on the occurrence
    # position of the guide: if upstream wrt indel position, guide position
    # correspond to pam start (easy adjustment), otherwise need to first remove
    # guide length (pam not considered)
    offset = guidepos if guidepos < guidepamlen else guidepos - (guidelen + pamlen - 1)
    # compute indel position within current guide candidate
    indelpos_rel = indelpos - offset
    if indelpos_rel < 0:  # overlapping insertion -> indel starts at position 0
        return 0
    return indelpos_rel


def filter_samples(
    samples: Set[str],
    variants: Dict[int, VariantRecord],
    position: int,
    indelpos: int,
    region_len: int,
    guidelen: int,
    pamlen: int,
    guidepos: int,
    i: int,
    chromcopy: int,
) -> Tuple[Set[str], int]:
    # compute indel position within current guide
    posrel = _compute_indel_posrel(indelpos, region_len, guidelen, pamlen, guidepos)
    # empty samples set or indel position -> do nothing
    if not samples or position == posrel:
        return samples, posrel
    try:
        samples_indel = variants[posrel].samples[i][chromcopy]  # samples carrying indel
    except KeyError:
        raise Exception
    return samples.intersection(samples_indel), posrel


def assign_na(values: List[str]) -> str:
    # return "NA" if the input values list is empty, otherwise construct a
    # comma-separated list
    if not values:
        return "NA"
    return ",".join(values)


def report_guide_haps(
    guides_samples: Dict[str, Set[str]],
    pos: int,
    positions: List[int],
    guides: List[str],
    samples: List[str],
    variants: List[str],
) -> Tuple[List[int], List[str], List[str], List[str]]:
    for guide in guides_samples:  # construct report lists
        positions.append(pos)  # guide position
        guides.append(guide)  # guide
        samples.append(
            assign_na(guides_samples[guide][0])
        )  # samples associated to the guide
        variants.append(
            assign_na(guides_samples[guide][1])
        )  # variants occurring in the guide
    return positions, guides, samples, variants


def find_snp_pos(
    guide: Guide, vmap: VariantMap
) -> Tuple[str, Dict[int, VariantRecord]]:
    # store variants found within the guide; look for positions within guide
    # sequence that carry variants
    variants = {
        i: vmap[guide.position + i]
        for i, _ in enumerate(guide)
        if guide.position + i in vmap
    }
    # compute guide reference sequence
    refseq = "".join(
        [variants[i].ref if i in variants else nt for i, nt in enumerate(guide)]
    )
    return refseq, variants


def retrieve_variant_samples(
    variants: Dict[int, VariantRecord], chromcopy: int
) -> Set[str]:
    return {
        sample
        for variant in variants.values()  # iterate over variants
        for samples_aa in variant.samples  # iterate over variant samples
        for sample in samples_aa[chromcopy]  # iterate over altallele samples
    }


def update_guideseq_snp(guideseq: List[str], position: int, alt: str) -> List[str]:
    # insert variant in guide candidate sequence (variants denoted by lowercase)
    return guideseq[:position] + [alt.lower()] + guideseq[(position + 1) :]


def update_sample_guide_snps(
    position: int, alt: str, guideseq: List[str], var_id: str, occurring_vars: List[str]
) -> Tuple[List[str], List[str]]:
    # update guide candidate sequence
    guideseq = update_guideseq_snp(guideseq, position, alt)
    occurring_vars.append(var_id)  # add current variant id
    return guideseq, occurring_vars


def filter_guides_snps(
    sample_guides_chromcopy: Dict[str, Tuple[List[str], List[str]]]
) -> Dict[str, Tuple[List[str], List[str]]]:
    # remove guides with no associated variant, and remove potential duplicate
    # ids from occurring variants list
    return {
        sample: (guide, list(dict.fromkeys(variants)))
        for sample, (guide, variants) in sample_guides_chromcopy.items()
        if variants
    }


def process_sample_guides_snps(
    refseq: str, samples: Set[str], variants: Dict[int, VariantRecord], chromcopy: int
) -> Dict[str, Tuple[List[str], List[str]]]:
    # initialize the sample-guide map with reference sequence for each sample
    # the ref guide is iteratively modified accordingly to the variants carried
    # by each sample
    # associate a list to store the variants contained in the reported guide
    # to each sample
    sample_guides = {sample: (list(refseq), []) for sample in samples}
    for pos, variant in variants.items():
        for i, samples_aa in enumerate(variant.samples):
            for sample in samples_aa[chromcopy]:
                guideseq, occurring_vars = sample_guides[sample]
                sample_guides[sample] = update_sample_guide_snps(
                    pos, variant.alt[i], guideseq, variant.id[i], occurring_vars
                )  # update sample's guide
    return sample_guides


def assign_sample_guide_snps(
    refseq: str, variants: Dict[int, VariantRecord], phased: bool
) -> Dict[int, Dict[str, Tuple[List[str], List[str]]]]:
    # if vcf is phased check variant occurrence on each chromosome copy
    chromcopies = 2 if phased else 1
    sample_guides = {i: {} for i in range(chromcopies)}
    for chromcopy in range(chromcopies):
        # retrieve samples carrying variants within input guide
        samples = retrieve_variant_samples(variants, chromcopy)
        if not samples:  # skip if no sample carries variant on current copy
            continue
        sample_guides[chromcopy] = filter_guides_snps(
            process_sample_guides_snps(refseq, samples, variants, chromcopy)
        )  # track guides with associated variants
    return sample_guides


def retrieve_guide_samples_snps(
    samples_guides: Dict[int, Dict[str, Tuple[List[str], List[str]]]],
    refseq: str,
):
    # initialize with reference guide -> no sample
    guide_dict = {refseq: (set(), [])}
    for samples_guides_map in samples_guides.values():  # reverse input map
        for sample, (guide, variants) in samples_guides_map.items():
            guide = "".join(guide)
            if guide not in guide_dict:
                guide_dict[guide] = (set(), variants)
            guide_dict[guide][0].add(sample)  # recover samples carrying this guide
    return guide_dict


def create_guide(
    guide: Guide, guideseq: str, samples: Set[str], variants: List[str]
) -> Guide:
    # hack to avoid side effects on the original guide's position (see init())
    position = guide.position if guide.right else guide.position + guide.guidelen
    g = Guide(
        position,
        list(guideseq),
        guide.guidelen,
        guide.pamlen,
        guide.strand,
        guide.debug,
        guide.right,
    )
    # set samples and variants for current guide
    if samples:
        g.set_samples(samples)
    if variants:
        g.set_variants(variants)
    return g


def reconstruct_haplotype_snps(guide: Guide, vmap: VariantMap) -> List[Guide]:
    # compute guide's reference sequence and positions carrying variants
    refseq, variants = find_snp_pos(guide, vmap)
    # map samples to their personal guide following haps if input vcf is phased
    samples_guides = assign_sample_guide_snps(refseq, variants, vmap.phased)
    # reverse the sample-guide dictionary to use guides as keys -> useful
    guides_samples = retrieve_guide_samples_snps(samples_guides, refseq)
    # update guide with samples and variants
    return [
        create_guide(guide, guideseq, *guides_samples[guideseq])
        for guideseq in guides_samples
    ]

def _isindel(variant: VariantRecord) -> bool:
    if variant.allelesnum == 1 and variant.vtype[0] == VTYPES[1]:
        return True  # indels are divided from potential other alleles at locus
    return False


def find_snp_indel_pos(
    region: IndelRegion, guide: Guide, vmap: VariantMap
) -> Tuple[str, Dict[int, VariantRecord]]:
    # store variants found within the guide; look for positions within guide
    # sequence that carry variants
    variants = {
        i: vmap[guide.position + i]
        for i, _ in enumerate(guide)
        if guide.position + i in vmap
    }
    # if guide starts from indel, adjust the guide's position
    if (0 in variants and _isindel(variants[0])) and region.indel_type == INDELTYPES[0]:
        guide.set_position(guide.position - 1)
    # retrieve guide reference
    start = guide.position
    stop = guide.position + len(guide) - region.indel_length if region.indel_type == INDELTYPES[0] else guide.position + len(guide) + region.indel_length
    refseq = "".join(region.sequence_ref[start:stop])
    return refseq, variants

def update_guideseq_indel(guideseq: List[str], position: int, ref: str, alt: str) -> List[str]:
    # insert variant in guide candidate sequence (variants denoted by lowercase)
    return guideseq[:position] + list(alt.lower()) + guideseq[position + len(ref):]


def update_sample_guide_indels(
    position: int, ref: str, alt: str, guideseq: List[str], var_id: str, occurring_vars: List[str],
) -> Tuple[List[str], List[str]]:
    # update guide candidate sequence
    guideseq = update_guideseq_indel(guideseq, position, ref, alt)
    occurring_vars.append(var_id)  # add current variant id
    return guideseq, occurring_vars


def process_sample_guides_indels(
    refseq: str, samples: Set[str], variants: Dict[int, VariantRecord], chromcopy: int
) -> Dict[str, Tuple[List[str], List[str]]]:
    # initialize the sample-guide map with reference sequence for each sample
    # the ref guide is iteratively modified accordingly to the variants carried
    # by each sample
    # associate a list to store the variants contained in the reported guide
    # to each sample
    sample_guides = {sample: (list(refseq), []) for sample in samples}
    variant_current = None
    for pos, variant in variants.items():
        if variant_current is not None and variant == variant_current:
            continue  # avoid repeated adjustments on insertions
        for i, samples_aa in enumerate(variant.samples):
            for sample in samples_aa[chromcopy]:
                guideseq, occurring_vars = sample_guides[sample]
                if pos < len(guideseq):
                    sample_guides[sample] = update_sample_guide_indels(
                        pos, variant.ref, variant.alt[i], guideseq, variant.id[i], occurring_vars
                    )  # update sample's guide
        variant_current = variant  # avoid repeated adjustments on insertions
    return sample_guides

def filter_guides_indels(
    sample_guides_chromcopy: Dict[str, Tuple[List[str], List[str]]], guidelen: int
) -> Dict[str, Tuple[List[str], List[str]]]:
    # remove guides with no associated variant, and remove potential duplicate
    # ids from occurring variants list
    return {
        sample: (guide, list(dict.fromkeys(variants)))
        for sample, (guide, variants) in sample_guides_chromcopy.items()
        if variants and len(guide) == guidelen
        # if (variants and indel.id[0] in variants) and len(guide) == guidelen
    }


def assign_sample_guide_indels(
    refseq: str, variants: Dict[int, VariantRecord], guidelen: int, phased: bool
) -> Dict[int, Dict[str, Tuple[List[str], List[str]]]]:
    # if vcf is phased check variant occurrence on each chromosome copy
    chromcopies = 2 if phased else 1
    sample_guides = {i: {} for i in range(chromcopies)}
    for chromcopy in range(chromcopies):
        # retrieve samples carrying variants within input guide
        samples = retrieve_variant_samples(variants, chromcopy)
        if not samples:  # skip if no sample carries variant on current copy
            continue
        sample_guides[chromcopy] = filter_guides_indels(
            process_sample_guides_indels(refseq, samples, variants, chromcopy), guidelen
        )  # track guides with associated variants
    return sample_guides



def reconstruct_haplotype_indels(region: IndelRegion, guide: Guide, vmap: VariantMap):
    # compute guide's reference sequence and positions carrying variants
    refseq, variants = find_snp_indel_pos(region, guide, vmap)
    # map samples to their personal guide following haps if input vcf is phased
    samples_guides = assign_sample_guide_indels(refseq, variants, len(guide), vmap.phased)
    # reverse the sample-guide dictionary to use guides as keys -> useful
    guides_samples = retrieve_guide_samples_snps(samples_guides, refseq)
    # update guide with samples and variants
    return [
        create_guide(guide, guideseq, *guides_samples[guideseq])
        for guideseq in guides_samples
        if len(guideseq) == len(guide)
    ]    

def track_haplotypes(
    region: Union[Region, IndelRegion], guides: List[Guide], vmap: VariantMap
) -> List[Guide]:
    if isinstance(region, IndelRegion):
        return [
            guide_
            for guide in guides
            for guide_ in reconstruct_haplotype_indels(region, guide, vmap)
        ]
    else:  # track haplotypes on SNP region
        return [
            guide_
            for guide in guides
            for guide_ in reconstruct_haplotype_snps(guide, vmap)
        ]
