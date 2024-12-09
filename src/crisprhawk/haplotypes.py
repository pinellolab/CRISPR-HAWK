"""
"""

from exception_handlers import exception_handler
from variants import VariantRecord, VTYPES

from typing import Tuple, List, Dict, Set

# iupac alphabet characters excluded canonical nucleotides [A, C, G, T, N]
IUPACSPECIALCHARS = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D"]


def find_variant_pos(guide: str, position: int, variants_maps: Dict[int, VariantRecord]) -> Tuple[str, Dict[int, VariantRecord]]:
    # store variants found within the guide; variants are denoted by iupac chars
    # different from the canonical nucleotides
    variants = {i: variants_maps[position + i] for i, nt in enumerate(guide) if nt in IUPACSPECIALCHARS}  
    # compute guide reference sequence
    refseq = "".join([variants[i].ref if i in variants else nt for i, nt in enumerate(guide)])
    return refseq, variants


def map_sample_to_guide(refseq: str, variants: Dict[int, VariantRecord], phased: bool) -> Dict[str, List[List[str]]]:
    # if vcf is phased check variant occurrence on each chromosome copy
    chromcopies = 2 if phased else 1  
    sample_guides = {i: {} for i in range(chromcopies)}
    for chromcopy in range(chromcopies): 
        # retrieve samples carrying variants within input guide
        samples = [s for v in variants.values() for aa in v.samples for s in aa[chromcopy]]
        if not samples:  # skip if no samples carries variant on current copy
            continue 
        # initialize the sample-guide map with reference sequence for each sample
        # the ref guide is iteratively modified accordingly to the variants carried 
        # by each sample
        sample_guides_chromcopy = {sample: list(refseq) for sample in samples}  
        for pos, variant in variants.items():  # iterate over variant positions 
            for i, aa in enumerate(variant.samples):  # iterate over alternative alleles
                if variant.vtype[i] == VTYPES[1]:  # indel -> skip
                    continue
                for sample in aa[chromcopy]:  # look for sample-specific variants
                    sample_seq = sample_guides_chromcopy[sample]
                    # insert variant in sample's guide - denote variants with
                    # lower case nucs
                    sample_seq[pos] = variant.alt[i].lower()  
                    sample_guides_chromcopy[sample] = sample_seq  # update sample's guide
        sample_guides[chromcopy] = sample_guides_chromcopy
    return sample_guides

def retrieve_guide_samples(sample_guides: Dict[str, List[List[str]]], refseq: str):
    guide_dict = {refseq: set()}  # initialize with reference guide -> no sample
    for chromcopy in sample_guides:  # reverse input map
        for sample, guide in sample_guides[chromcopy].items():  
            guide = "".join(guide)
            if guide not in guide_dict:
                guide_dict[guide] = set()
            guide_dict[guide].add(sample)  # recover samples carrying this guide
    return guide_dict

def report_guide_haps(guides_samples: Dict[str, Set[str]], pos: int, positions: List[int], guides: List[str], samples: List[str]) -> Tuple[List[int], List[str], List[str]]:
    for guide in guides_samples:   # construct report lists
        positions.append(pos)
        guides.append(guide)
        samples.append(",".join(sorted(guides_samples[guide])))  # sort for reproducibility
    return positions, guides, samples


def reconstruct_guide_haps(guides: List[str], positions: List[int], strand: int, guidelen: int, variants_maps: Dict[int, VariantRecord], phased: bool) -> Tuple[List[int], List[str], List[str]]:
    # report fields: position, guide, samples
    positions_rep, guides_rep, samples_rep = [], [], []
    # reconstruct haps for each guide and retrieve associated samples
    for i, guide in enumerate(guides):  
        # recover guide start position
        pos = positions[i] - guidelen if strand == 0 else positions[i] 
        # compute reference guide sequence and positions with variants
        refseq, variants = find_variant_pos(guide, pos, variants_maps)
        # map samples to their personal guide following haps if input vcf is phased
        samples_guides = map_sample_to_guide(refseq, variants, phased)
        # reverse the sample-guide dictionary to use guides as keys -> useful
        # for reporting purposes
        guides_samples = retrieve_guide_samples(samples_guides, refseq)
        # construct report for current strand
        positions_rep, guides_rep, samples_rep = report_guide_haps(guides_samples, positions[i], positions_rep, guides_rep, samples_rep)
    return positions_rep, guides_rep, samples_rep

def track_haplotypes(guides: Tuple[Tuple[List[int], List[int]], Tuple[List[str], List[str]]], variants_maps: Dict[int, VariantRecord], guidelen: int, phased: bool) -> Tuple[Tuple[List[int], List[int]], Tuple[List[str], List[str]], Tuple[List[str], List[str]]]:
    # track haplotypes on forward and reverse strands
    # if phased, check variants copy mapping
    position_proc, guides_seq_proc, samples_proc = [[], []], [[], []], [[], []]
    for strand in [0, 1]:  # 0 is + strand, 1 is - strand
        # reconstruct guides haplotypes 
        position_proc[strand], guides_seq_proc[strand], samples_proc[strand] = reconstruct_guide_haps(guides[1][strand], guides[0][strand], strand, guidelen, variants_maps, phased)
    return (position_proc[0], position_proc[1]), (guides_seq_proc[0], guides_seq_proc[1]), (samples_proc[0], samples_proc[1])
