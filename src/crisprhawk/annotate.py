"""
"""

from exception_handlers import exception_handler
from scores import azimuth
from guide import Guide, GUIDESEQPAD

from hapsolver import Region
from typing import List, Dict

import numpy as np

import pickle
import os

def reverse_guides(guides: List[Guide]) -> List[Guide]:
    # compute reverse complement sequence for guides occurring on reverse strand
    for guide in guides:
        if guide.strand == 1:  # guide on reverse strand
            guide.reverse_complement()
    return guides

def azimuth_score(guides: List[Guide], debug: bool) -> List[Guide]:
    # create guides np.ndarray required by azimuth; each guide must have 4 nts
    # upstream the guide sequence, and 3 nts downstream the pam
    guides_seqs = np.array([guide.sequence[(GUIDESEQPAD - 4):(-GUIDESEQPAD + 3)].upper() for guide in guides])
    try:  # compute azimuth scores
        azimuth_scores = azimuth(guides_seqs)  
    except Exception as e:
        exception_handler(ValueError, "Azimuth score calculation failed", os.EX_DATAERR, debug, e)
    assert len(azimuth_scores) == len(guides)  # should match
    for i, score in enumerate(azimuth_scores):
        guides[i].set_azimuth_score(score)  # assign score to each guide
    return guides

def compute_position(guides: List[Guide], region_start: int) -> List[Guide]:
    for guide in guides:  # recover genomic position of guides
        position = guide.position + region_start + 1  # from 0-based to 1-based
        guide.set_position(position)
    return guides

def annotate_variants(guides: List[Guide], debug: bool) -> List[Guide]:
    for guide in guides:
        guidepamlen = guide.guidelen + guide.pamlen  # total guide length
        guide_vars = set()  # variants occurring in variant
        for variant in guide.variants.split(","):
            try:  # retrieve each variant position
                variant_position = int(variant.split("-")[1])
            except TypeError as e:
                exception_handler(TypeError, f"Variant {variant} seems to have a non int position", os.EX_DATAERR, debug, e)
            # assess whether the snp occurs within the guide or is part of the haplotype
            if guide.position <= variant_position <= guide.position + guidepamlen:
                guide_vars.add(variant)
        # set variants ids to current guide
        guide_vars = ",".join(sorted(guide_vars)) if guide_vars else "NA"
        guide.set_variants(guide_vars)
    return guides


def annotate_guides(guides: Dict[Region, List[Guide]], debug: bool):
    # annotate guides with scores, variants and adjust positions
    for region, guides_list in guides.items():
        # compute reverse complement for guides occurring on rev strand
        guides_list = reverse_guides(guides_list)
        # compute guides position in genome
        guides_list = compute_position(guides_list, region.start)
        # set variants for current guide
        guides_list = annotate_variants(guides_list, debug)
        # annotate each guide with azimuth scores
        guides[region] = azimuth_score(guides_list, debug) 
    return guides        
        

    

    