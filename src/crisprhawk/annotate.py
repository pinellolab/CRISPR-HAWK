""" """

from exception_handlers import exception_handler
from scores import azimuth
from guide import Guide, GUIDESEQPAD
from utils import print_verbosity, VERBOSITYLVL

from hapsolver import Region
from typing import List, Dict
from time import time

import numpy as np

import pickle
import os


def reverse_guides(guides: List[Guide], verbosity: int) -> List[Guide]:
    # compute reverse complement sequence for guides occurring on reverse strand
    print_verbosity(
        "Reversing guides occurring on reverse strand", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # reversal start time
    for guide in guides:
        if guide.strand == 1:  # guide on reverse strand
            guide.reverse_complement()
    print_verbosity(
        f"Guides reversed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def azimuth_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    # create guides np.ndarray required by azimuth; each guide must have 4 nts
    # upstream the guide sequence, and 3 nts downstream the pam
    print_verbosity("Computing Azimuth score", verbosity, VERBOSITYLVL[3])
    start = time()  # azimuth score start time
    guides_seqs = np.array(
        [
            guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
            for guide in guides
        ]
    )
    try:  # compute azimuth scores
        azimuth_scores = azimuth(guides_seqs)
    except Exception as e:
        exception_handler(
            ValueError, "Azimuth score calculation failed", os.EX_DATAERR, debug, e
        )
    assert len(azimuth_scores) == len(guides)  # should match
    for i, score in enumerate(azimuth_scores):
        guides[i].set_azimuth_score(score)  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def compute_position(
    guides: List[Guide], region_start: int, verbosity: int
) -> List[Guide]:
    print_verbosity("Computing guides genomic positions", verbosity, VERBOSITYLVL[3])
    start = time()  # position calculation start time
    for guide in guides:  # recover genomic position of guides
        position = guide.position + region_start + 1  # from 0-based to 1-based
        guide.set_position(position)
    print_verbosity(
        f"Genomic positions computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def annotate_variants(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity(
        "Annotating variants occurring in guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # position calculation start time
    for guide in guides:
        guidepamlen = guide.guidelen + guide.pamlen  # total guide length
        guide_vars = set()  # variants occurring in variant
        for variant in guide.variants.split(","):
            try:  # retrieve each variant position
                variant_position = int(variant.split("-")[1])
            except TypeError as e:
                exception_handler(
                    TypeError,
                    f"Variant {variant} seems to have a non int position",
                    os.EX_DATAERR,
                    debug,
                    e,
                )
            # assess whether the snp occurs within the guide or is part of the haplotype
            if guide.position <= variant_position <= guide.position + guidepamlen:
                guide_vars.add(variant)
        # set variants ids to current guide
        guide_vars = ",".join(sorted(guide_vars)) if guide_vars else "NA"
        guide.set_variants(guide_vars)
    print_verbosity(
        f"Variants annotated in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def annotate_guides(
    guides: Dict[Region, List[Guide]], verbosity: int, debug: bool
) -> Dict[Region, List[Guide]]:
    # annotate guides with scores, variants and adjust positions
    print_verbosity("Annotating guides", verbosity, VERBOSITYLVL[1])
    start = time()  # annotation start time
    for region, guides_list in guides.items():
        # compute reverse complement for guides occurring on rev strand
        guides_list = reverse_guides(guides_list, verbosity)
        # compute guides position in genome
        guides_list = compute_position(guides_list, region.start, verbosity)
        # set variants for current guide
        guides_list = annotate_variants(guides_list, verbosity, debug)
        # annotate each guide with azimuth scores
        guides[region] = azimuth_score(guides_list, verbosity, debug)
    print_verbosity(
        f"Annotation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides
