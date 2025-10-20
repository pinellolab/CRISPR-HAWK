"""
This module provides functions and utilities for annotating CRISPR guide RNAs with
various scores, sequence features, variant information, and functional or gene
annotations.

It supports the calculation of efficiency and specificity scores, variant annotation,
allele frequency assignment, GC content, out-of-frame scores, and integration with
external tools for off-target prediction. The module is designed to process and
enrich guide RNA data for downstream analysis in genome editing applications.
"""

from .crisprhawk_error import (
    CrisprHawkCfdScoreError,
    CrisprHawkAzimuthScoreError,
    CrisprHawkRs3ScoreError,
    CrisprHawkDeepCpf1ScoreError,
    CrisprHawkAnnotationError,
    CrisprHawkGcContentError,
    CrisprHawkOOFrameScoreError,
)
from .config_crispritz import CrispritzConfig
from .exception_handlers import exception_handler
from .scores import azimuth, cfdon, rs3, deepcpf1, elevationon, ooframe_score
from .bedfile import BedAnnotation
from .guide import Guide, GUIDESEQPAD
from .utils import print_verbosity, flatten_list, VERBOSITYLVL
from .region import Region
from .pam import PAM, CPF1, SPCAS9, XCAS9
from .offtargets import search_offtargets

from collections import defaultdict
from typing import List, Dict, Union, Set, Tuple
from Bio.SeqUtils import gc_fraction
from time import time

import numpy as np

import os

ANNDIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "annotations"
)  # annotation data directory


def reverse_guides(guides: List[Guide], verbosity: int) -> List[Guide]:
    """Reverses the sequence of guides that are located on the reverse strand.

    This function computes the reverse complement for each guide on the reverse
    strand and updates their sequences accordingly.

    Args:
        guides (List[Guide]): List of Guide objects to process.
        verbosity (int): Verbosity level for logging.

    Returns:
        List[Guide]: The list of guides with reverse strand guides reversed.
    """
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
    """Computes Azimuth scores for a list of guide RNAs.

    This function calculates the Azimuth efficiency score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with Azimuth scores assigned.
    """
    # create guides np.ndarray required by azimuth; each guide must have 4 nts
    # upstream the guide sequence, and 3 nts downstream the pam
    if not guides:
        return guides
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
            CrisprHawkAzimuthScoreError,
            "Azimuth score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(azimuth_scores) == len(guides)  # should match
    for i, score in enumerate(azimuth_scores):
        guides[i].azimuth_score = score  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def rs3_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes RS3 scores for a list of guide RNAs.

    This function calculates the RS3 efficiency score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with RS3 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing RS3 score", verbosity, VERBOSITYLVL[3])
    start = time()  # rs3 score start time
    guides_seqs = [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]
    try:  # compute rs3 scores
        rs3_scores = rs3(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkRs3ScoreError,
            "RS3 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(rs3_scores) == len(guides)  # should match
    for i, score in enumerate(rs3_scores):
        guides[i].rs3_score = score  # assign score to each guide
    print_verbosity(
        f"RS3 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def deepcpf1_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes DeepCpf1 scores for a list of guide RNAs.

    This function calculates the DeepCpf1 efficiency score for each guide and
    updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with DeepCpf1 scores assigned.
    """
    if not guides:
        return guides
    print_verbosity("Computing DeepCpf1 score", verbosity, VERBOSITYLVL[3])
    start = time()  # deepcpf1 score start time
    guides_seqs = [
        guide.sequence[(GUIDESEQPAD - 4) : (-GUIDESEQPAD + 3)].upper()
        for guide in guides
    ]
    try:  # compute deepcpf1 scores
        deepcpf1_scores = deepcpf1(guides_seqs)
    except Exception as e:
        exception_handler(
            CrisprHawkDeepCpf1ScoreError,
            "DeepCpf1 score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    assert len(deepcpf1_scores) == len(guides)  # should match
    for i, score in enumerate(deepcpf1_scores):
        guides[i].deepcpf1_score = score  # assign score to each guide
    print_verbosity(
        f"DeepCpf1 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def group_guides_position(
    guides: List[Guide], debug: bool
) -> Dict[str, Tuple[Union[None, Guide], List[Guide]]]:
    """Groups guides by their genomic position and strand.

    This function organizes guides into groups based on their start position and
    strand, identifying the reference guide and associated variant guides for each
    group.

    Args:
        guides (List[Guide]): List of Guide objects to group.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Dict[str, Tuple[Union[None, Guide], List[Guide]]]: A dictionary mapping
            position keys to tuples containing the reference guide and a list of
            guides at that position.
    """

    class _GuideGroup:
        """Helper class to group guides by position and strand.

        This class stores a reference guide and a list of guides for a specific
        genomic position and strand.
        """

        def __init__(self) -> None:
            self._refguide = None
            self._guides = []

        def to_tuple(self) -> Tuple[Union[Guide, None], List[Guide]]:
            return (self._refguide, self._guides)

    pos_guide = defaultdict(_GuideGroup)
    for guide in guides:
        poskey = f"{guide.start}_{guide.strand}"
        if guide.samples == "REF":  # reference guide
            if pos_guide[poskey]._refguide is not None:
                exception_handler(
                    CrisprHawkCfdScoreError,
                    f"Duplicate REF guide at position {guide.start}? CFDon/Elevation-on calculation failed",
                    os.EX_DATAERR,
                    debug,
                )
            pos_guide[poskey]._refguide = guide  # type: ignore
        pos_guide[poskey]._guides.append(guide)
    return {poskey: g.to_tuple() for poskey, g in pos_guide.items()}


def cfdon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes CFDon scores for a list of guide RNAs.

    This function calculates the CFDon specificity score for each guide and updates
    the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with CFDon scores assigned.
    """
    print_verbosity("Computing CFDon score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    for _, (guide_ref, guides_g) in guide_groups.items():
        try:
            cfdon_scores = cfdon(guide_ref, guides_g, debug)
        except Exception as e:
            exception_handler(
                CrisprHawkCfdScoreError,
                "CFDon score calculation failed",
                os.EX_DATAERR,
                debug,
                e,
            )
        for i, score in enumerate(cfdon_scores):
            guides_g[i].cfdon_score = score
    # revert grouped guides by position into list
    guides = flatten_list([guides_g for _, (_, guides_g) in guide_groups.items()])
    print_verbosity(
        f"CFDon scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def elevationon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes Elevation-on scores for a list of guide RNAs.

    This function calculates the Elevation-on efficiency score for each guide and
    updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to score.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with Elevation-on scores assigned.
    """
    print_verbosity("Computing Elevation-on score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    guides = elevationon(guide_groups)
    print_verbosity(
        f"Elevation-on scores computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def polish_variants_annotation(guide: Guide, variants: Set[str]) -> Set[str]:
    """Validates and filters variants that overlap with a guide sequence.

    This function checks each variant to ensure it matches the expected sequence
    context within the guide and returns only validated variants.

    Args:
        guide (Guide): The guide object whose sequence is used for validation.
        variants (Set[str]): Set of variant identifiers to validate.

    Returns:
        Set[str]: The set of validated variant identifiers that match the guide
            sequence.
    """
    # map variant positions to variant strings for quick lookup
    varposmap = {int(variant.split("-")[1]): variant for variant in variants}
    validated_variants = set()
    # check each nucleotide position in the guide sequence
    for seqidx, nt in enumerate(guide.guidepam):
        pos = guide.start + seqidx
        if pos not in varposmap:
            continue
        # parse variant information: "chrom-pos-ref/alt"
        variant_id = varposmap[pos]
        ref, alt = variant_id.split("-")[2].split("/")
        if len(ref) != len(
            alt
        ):  # patch to handle indels overlapping guide's start or stop
            validated_variants.add(variant_id)
            continue
        # calculate sequence length difference (indel offset)
        offset = max(0, len(alt) - len(ref))
        # validate: lowercase sequence in guide should match uppercase alt allele
        guide_segment = guide.guidepam[seqidx : seqidx + offset + 1]
        if guide_segment.islower() and guide_segment.upper() == alt:
            validated_variants.add(variant_id)
    return validated_variants


def annotate_variants(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Annotates guides with variants that overlap their sequence.

    This function identifies and validates variants that occur within each guide's
    sequence and updates the guide objects with the relevant variant information.

    Args:
        guides (List[Guide]): List of Guide objects to annotate.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with annotated variant information.
    """
    guides_lst = []  # reported guides
    print_verbosity(
        "Annotating variants occurring in guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # position calculation start time
    for guide in guides:
        guide_vars = set()  # variants occurring in variant
        is_reference = False
        for variant in guide.variants.split(","):
            if variant == "NA":  # reference sequence (no variants)
                if guide_vars:  # shouldn't have both NA and actual variants
                    exception_handler(
                        ValueError, "Forbidden NA variant", os.EX_DATAERR, debug
                    )
                guide_vars.add("NA")
                is_reference = True
                break
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
            if guide.start <= variant_position < guide.stop:
                guide_vars.add(variant)
        if guide_vars:  # process guides with variants
            if not is_reference:  # polish variants to validate sequence matches
                guide_vars = polish_variants_annotation(guide, guide_vars)
            guide.variants = ",".join(sorted(guide_vars))
            guides_lst.append(guide)
    print_verbosity(
        f"Variants annotated in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides_lst


def annotate_variants_afs(guides: List[Guide], verbosity: int) -> List[Guide]:
    """Annotates guides with allele frequencies for variants in their sequence.

    This function assigns allele frequency values to each guide based on the
    variants present in its sequence.

    Args:
        guides (List[Guide]): List of Guide objects to annotate.
        verbosity (int): Verbosity level for logging.

    Returns:
        List[Guide]: The list of guides with annotated allele frequency information.
    """
    guides_lst = []  # reported guides
    print_verbosity(
        "Annotating variants allele frequencies in guides", verbosity, VERBOSITYLVL[3]
    )
    start = time()  # position calculation start time
    for guide in guides:
        afs = (
            [
                str(guide.afs[v]) if str(guide.afs[v]) != "nan" else "NA"
                for v in guide.variants.split(",")
            ]
            if guide.variants != "NA"
            else ["NA"]
        )
        guide.afs_str = afs
        guides_lst.append(guide)
    print_verbosity(
        f"Variants allele frequencies annotated in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides_lst


def _funcann(
    guide: Guide, bedannotation: BedAnnotation, contig: str, debug: bool
) -> Guide:
    """Annotates a guide with functional features from a BED annotation.

    This function fetches annotation features overlapping the guide and assigns
    the relevant annotation to the guide object.

    Args:
        guide (Guide): The guide object to annotate.
        bedannotation (BedAnnotation): The BED annotation object to query.
        contig (str): The contig or chromosome name.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Guide: The guide object with functional annotation set.
    """
    try:  # fetch annotation features overlapping input guide
        annotation = bedannotation.fetch_features(contig, guide.start, guide.stop)
    except Exception as e:
        exception_handler(
            CrisprHawkAnnotationError,
            f"Guides annotation failed on {guide}",
            os.EX_DATAERR,
            debug,
            e,
        )
    # if no annotation, return NA value; annotation values on 4th BED column
    annotation = ",".join([e.split()[3] for e in annotation]) if annotation else "NA"
    guide.funcann = annotation
    return guide


def _retrieve_gene_name(field: str) -> str:
    """Extracts the gene name from a semicolon-separated annotation field.

    This function searches for the 'gene_name=' substring and returns the corresponding
    gene name if present.

    Args:
        field (str): The annotation field string to search.

    Returns:
        str: The extracted gene name, or an empty string if not found.
    """
    i = field.find("gene_name=")
    return "" if i == -1 else field[i + 10 : field.find(";", i + 10)]


def _geneann(
    guide: Guide, bedannotation: BedAnnotation, contig: str, debug: bool
) -> Guide:
    """Annotates a guide with gene features from a BED annotation.

    This function fetches gene annotation features overlapping the guide and
    assigns the relevant gene annotation to the guide object.

    Args:
        guide (Guide): The guide object to annotate.
        bedannotation (BedAnnotation): The BED annotation object to query.
        contig (str): The contig or chromosome name.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Guide: The guide object with gene annotation set.
    """
    try:  # fetch annotation features overlapping input guide
        annotation = bedannotation.fetch_features(contig, guide.start, guide.stop)
    except Exception as e:
        exception_handler(
            CrisprHawkAnnotationError,
            f"Guides gene annotation failed on {guide}",
            os.EX_DATAERR,
            debug,
            e,
        )
    # if no annotation, return NA value; annotation values on 4th BED column
    annotation = (
        ",".join(
            [
                f"{fields[7]}:{_retrieve_gene_name(fields[9])}"
                for e in annotation
                for fields in [e.split()]
            ]
        )
        if annotation
        else "NA"
    )
    guide.geneann = annotation
    return guide


def ann_guides(
    guides: List[Guide],
    contig: str,
    annotations: List[str],
    atype: int,  # 0 -> regular annotation; 1 -> gene annotation
    verbosity: int,
    debug: bool,
) -> List[Guide]:
    """Annotates guides with functional or gene features from BED files.

    This function applies either regular or gene annotation to each guide using
    the provided BED annotation files.

    Args:
        guides (List[Guide]): List of Guide objects to annotate.
        contig (str): The contig or chromosome name.
        annotations (List[str]): List of BED annotation file paths.
        atype (int): Annotation type (0 for regular annotation, 1 for gene annotation).
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with applied annotations.
    """
    print_verbosity("Starting guides annotation", verbosity, VERBOSITYLVL[3])
    start = time()  # functional annotation start time
    assert atype in {0, 1}  # used to set the proper field in guides
    # idx = 22 if atype == 1 else 3
    guides_ann = []
    for fann in annotations:
        bedann = BedAnnotation(fann, verbosity, debug)  # load annotation bed
        guides_ann = [
            (
                _funcann(guide, bedann, contig, debug)
                if atype == 0
                else _geneann(guide, bedann, contig, debug)
            )
            for guide in guides
        ]
    assert len(guides) == len(guides_ann)  # type: ignore
    print_verbosity(
        f"Guides functional annotation completed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides_ann or guides


def gc_content(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Computes the GC content for each guide RNA sequence.

    This function calculates the GC content (excluding the PAM) for each guide
    and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to process.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with GC content assigned.
    """
    print_verbosity("Computing GC content", verbosity, VERBOSITYLVL[3])
    start = time()  # GC content calculation start time
    try:  # compute gc content (PAM excluded)
        for guide in guides:
            guide.gc = gc_fraction(guide.guide)
    except Exception as e:
        exception_handler(
            CrisprHawkGcContentError,
            "GC content calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    print_verbosity(
        f"GC content computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def outofframe_score(
    guides: List[Guide], guidelen: int, right: bool, verbosity: int, debug: bool
) -> List[Guide]:
    """Computes the out-of-frame score for each guide RNA sequence.

    This function calculates the likelihood that a guide induces an out-of-frame
    mutation and updates the guide objects with the computed values.

    Args:
        guides (List[Guide]): List of Guide objects to process.
        guidelen (int): Length of the guide sequence.
        right (bool): Whether the guide is extracted downstream (right side) of
            the PAM.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        List[Guide]: The list of guides with out-of-frame scores assigned.
    """
    print_verbosity("Computing out-of-frame score", verbosity, VERBOSITYLVL[3])
    start = time()  # out-of-frame score calculation start time
    try:  # compute out-of-frame score
        idx = GUIDESEQPAD if right else GUIDESEQPAD + guidelen
        # scores = [0] * len(guides)
        scores = ooframe_score(guides, idx)
    except Exception as e:
        exception_handler(
            CrisprHawkOOFrameScoreError,
            "Out-of-frame score calculation failed",
            os.EX_DATAERR,
            debug,
            e,
        )
    for i, score in enumerate(scores):  # set out-of-frame score for each guide
        guides[i].ooframe_score = score
    print_verbosity(
        f"Out-of-frame score computed in {time() - start:.2f}s",
        verbosity,
        VERBOSITYLVL[3],
    )
    return guides


def annotate_guides(
    guides: Dict[Region, List[Guide]],
    annotations: List[str],
    gene_annotations: List[str],
    pam: PAM,
    compute_elevation: bool,
    estimate_offtargets: bool,
    crispritz_config: Union[None, CrispritzConfig],
    mm: int,
    bdna: int,
    brna: int,
    offtargets_annotations: List[str],
    offtargets_annotation_colnames: List[str],
    crispritz_index: str,
    guidelen: int,
    right: bool,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    """Annotates guides with scores, sequence features, and functional or gene
    information.

    This function processes each region's guides by annotating them with variant,
    efficiency, specificity, and functional or gene features, as well as off-target
    predictions if requested.

    Args:
        guides (Dict[Region, List[Guide]]): Dictionary mapping regions to lists
            of Guide objects.
        annotations (List[str]): List of BED annotation file paths for functional
            annotation.
        gene_annotations (List[str]): List of BED annotation file paths for gene
            annotation.
        pam (PAM): PAM object specifying the CRISPR system.
        compute_elevation (bool): Whether to compute Elevation-on scores.
        estimate_offtargets (bool): Whether to estimate off-targets using CRISPRitz.
        crispritz_config (Union[None, CrispritzConfig]): CRISPRitz configuration
            object.
        mm (int): Maximum number of mismatches for off-target search.
        bdna (int): Maximum number of DNA bulges for off-target search.
        brna (int): Maximum number of RNA bulges for off-target search.
        crispritz_index (str): Path to the CRISPRitz genome index.
        guidelen (int): Length of the guide sequence.
        right (bool): Whether the guide is extracted downstream (right side) of
            the PAM.
        outdir (str): Output directory for results.
        threads (int): Number of threads to use.
        verbosity (int): Verbosity level for logging.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Dict[Region, List[Guide]]: The dictionary of regions with annotated guides.
    """
    # annotate guides with scores, variants and adjust positions
    print_verbosity("Annotating guides", verbosity, VERBOSITYLVL[1])
    start = time()  # annotation start time
    for region, guides_list in guides.items():
        # set variants for current guide
        guides_list = annotate_variants(guides_list, verbosity, debug)
        # add allele frequencies for variants occurring in guides
        guides_list = annotate_variants_afs(guides_list, verbosity)
        # compute reverse complement for guides occurring on rev strand
        guides_list = reverse_guides(guides_list, verbosity)
        if pam.cas_system in [SPCAS9, XCAS9]:  # cas9 system pam
            # annotate each guide with azimuth scores
            guides_list = azimuth_score(guides_list, verbosity, debug)
            # annotate each guide with rs3 scores
            guides_list = rs3_score(guides_list, verbosity, debug)
            # annotate each guide with CFDon scores
            guides_list = cfdon_score(guides_list, verbosity, debug)
        if pam.cas_system == CPF1:  # cpf1 system pam
            guides_list = deepcpf1_score(guides_list, verbosity, debug)
        if compute_elevation and (guidelen + len(pam) == 23 and not right):
            # elevation requires 23 bp long sequences, where last 3 bp are pam
            guides_list = elevationon_score(guides_list, verbosity, debug)
        # compute gc content (pam excluded) for each guide
        guides_list = gc_content(guides_list, verbosity, debug)
        # compute out-of-frame score
        guides_list = outofframe_score(guides_list, guidelen, right, verbosity, debug)
        if annotations:  # annotate each guide
            guides_list = ann_guides(
                guides_list,
                region.contig,
                annotations,
                0,
                verbosity,
                debug,
            )
        if gene_annotations:  # annotate each guide with gene data
            guides_list = ann_guides(
                guides_list, region.contig, gene_annotations, 1, verbosity, debug
            )
        if estimate_offtargets:  # estimate off-targets for each guide
            assert crispritz_config  # shouldn't be None if we got here
            guides_list = search_offtargets(
                guides_list,
                pam,
                crispritz_index,
                region,
                crispritz_config,
                mm,
                bdna,
                brna,
                offtargets_annotations,
                offtargets_annotation_colnames,
                guidelen,
                compute_elevation,
                right,
                threads,
                outdir,
                verbosity,
                debug,
            )
        guides[region] = guides_list  # store annotated guides
    print_verbosity(
        f"Annotation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides
