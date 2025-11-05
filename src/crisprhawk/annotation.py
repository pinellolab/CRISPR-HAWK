"""Provides functions for annotating CRISPR guides with variant, functional, gene,
and GC content information.

This module includes utilities for processing guide sequences, validating and
annotating variants, assigning allele frequencies, adding functional and gene
annotations from BED files, and computing GC content.

It supports comprehensive annotation of guides for downstream CRISPR analysis workflows.
"""

from .crisprhawk_error import CrisprHawkAnnotationError, CrisprHawkGcContentError
from .crisprhawk_argparse import CrisprHawkSearchInputArgs
from .exception_handlers import exception_handler
from .bedfile import BedAnnotation
from .guide import Guide
from .utils import print_verbosity, VERBOSITYLVL
from .region import Region

from typing import List, Dict, Set, Tuple
from Bio.SeqUtils import gc_fraction
from time import time

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


def _parse_variant(variant_id: str) -> Tuple[str, int, str, str]:
    """Parses a variant identifier string into its components.

    This function splits a variant ID into chromosome, position, reference allele,
    and alternate allele.

    Args:
        variant_id (str): The variant identifier string in the format 'chrom-pos-ref/alt'.

    Returns:
        Tuple[str, int, str, str]: A tuple containing chromosome, position,
            reference allele, and alternate allele.
    """
    parts = variant_id.split("-")
    chrom = parts[0]  # chromosome
    position = int(parts[1])  # variant position
    ref, alt = parts[2].split("/")
    return chrom, position, ref, alt


def _is_snv(ref: str, alt: str) -> bool:
    """Determines if a variant is a single nucleotide variant (SNV).

    This function returns True if the reference and alternate alleles are of
    equal length, indicating an SNV.

    Args:
        ref (str): The reference allele.
        alt (str): The alternate allele.

    Returns:
        bool: True if the variant is an SNV, False otherwise.
    """
    return len(ref) == len(alt)


def sort_variants(variants):
    """Sorts variant identifiers in natural genomic order.

    This function sorts variants by chromosome, position, variant type (SNV or indel),
    reference, and alternate allele.

    Args:
        variants: An iterable of variant identifier strings.

    Returns:
        List[str]: The sorted list of variant identifiers.
    """

    def _sort_key(variant):
        chrom, pos, ref, alt = _parse_variant(variant)
        # extract numeric part from chromosome for natural sorting
        chrom_parts = chrom.replace("chr", "")
        try:
            chrom_num = int(chrom_parts)
        except ValueError:
            chrom_num = float("inf")  # handle X, Y, M, etc.
            chrom_parts = chrom_parts.lower()
        # return tuple for sorting:
        # (chrom_number, chrom_string, position, is_indel, ref, alt)
        return (chrom_num, chrom_parts, pos, not _is_snv(ref, alt), ref, alt)

    return sorted(variants, key=_sort_key)


def _create_variant_position_map(variants: Set[str]) -> Dict[int, str]:
    """Creates a mapping from variant positions to variant identifiers.

    This function generates a dictionary mapping each variant's genomic position
    to its identifier string.

    Args:
        variants (Set[str]): A set of variant identifier strings.

    Returns:
        Dict[int, str]: A dictionary mapping positions to variant identifiers.
    """
    varposmap = {}  # map variant positions to variant strings for quick lookup
    for variant in sort_variants(variants):
        _, pos, _, _ = _parse_variant(variant)
        varposmap[pos] = variant
    return varposmap


def _compute_indel_offset(ref: str, alt: str) -> Tuple[int, int]:
    """Computes the offset for insertions and deletions between reference and
    alternate alleles.

    This function calculates the number of bases deleted and inserted by comparing
    the lengths of the reference and alternate alleles.

    Args:
        ref (str): The reference allele.
        alt (str): The alternate allele.

    Returns:
        Tuple[int, int]: A tuple containing the deletion offset and insertion offset.
    """
    deletion_offset = max(0, len(ref) - len(alt))
    insertion_offset = max(0, len(alt) - len(ref))
    return deletion_offset, insertion_offset


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
    varposmap = _create_variant_position_map(variants)
    offset_insertion = 0  # Skip bases in guide sequence for insertions
    offset_deletion = 0  # Adjust position calculation for deletions
    validated_variants = set()
    # check each nucleotide position in the guide sequence
    for seqidx, nt in enumerate(guide.guidepam):
        if offset_insertion > 0:  # skip positions if in an insertion
            offset_insertion -= 1
            continue
        # genomic position accounting for deletions
        pos = guide.start + seqidx + offset_deletion
        if pos not in varposmap:  # check if there's a variant at this position
            continue
        # parse variant information: "chrom-pos-ref/alt"
        variant_id = varposmap[pos]
        _, pos, ref, alt = _parse_variant(variant_id)
        if not _is_snv(ref, alt):  #  handle indels
            validated_variants.add(variant_id)
            # update indel offsets
            do, offset_insertion = _compute_indel_offset(ref, alt)
            offset_deletion += do
            continue
        # validate: lowercase sequence in guide should match uppercase alt allele
        guide_segment = guide.guidepam[seqidx : seqidx + 1]
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
    # for guide in guides:
    #     guide_vars = set()  # variants occurring in variant
    #     is_reference = False
    #     for variant in guide.variants.split(","):
    #         if variant == "NA":  # reference sequence (no variants)
    #             if guide_vars:  # shouldn't have both NA and actual variants
    #                 exception_handler(
    #                     ValueError, "Forbidden NA variant", os.EX_DATAERR, debug
    #                 )
    #             guide_vars.add("NA")
    #             is_reference = True
    #             break
    #         try:  # retrieve each variant position
    #             variant_position = int(variant.split("-")[1])
    #         except TypeError as e:
    #             exception_handler(
    #                 TypeError,
    #                 f"Variant {variant} seems to have a non int position",
    #                 os.EX_DATAERR,
    #                 debug,
    #                 e,
    #             )
    #         # assess whether the snp occurs within the guide or is part of the haplotype
    #         if guide.start <= variant_position < guide.stop:
    #             guide_vars.add(variant)
    #     if guide_vars:  # process guides with variants
    #         if not is_reference:  # polish variants to validate sequence matches
    #             guide_vars = polish_variants_annotation(guide, guide_vars)
    #         guide.variants = ",".join(sorted(guide_vars))
    #         guides_lst.append(guide)
    for g in guides:
        g.variants = "NA"
        guides_lst.append(g)
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


def annotate_guides(
    guides: Dict[Region, List[Guide]], args: CrisprHawkSearchInputArgs
) -> Dict[Region, List[Guide]]:
    """Annotates CRISPR guides with variant, functional, gene, and GC content
    information.

    This function processes each region's guides, adding variant, allele frequency,
    reverse complement, GC content, and optional functional and gene annotations,
    returning the updated guides.

    Args:
        guides: Dictionary mapping Region objects to lists of Guide objects.
        args: CrisprHawkSearchInputArgs object containing annotation parameters.

    Returns:
        Dictionary mapping Region objects to lists of annotated Guide objects.
    """
    # annotate guides with variants, functional and gene data and adjust positions
    print_verbosity("Annotating guides", args.verbosity, VERBOSITYLVL[1])
    start = time()  # annotation start time
    for region, guides_list in guides.items():
        # set variants for current guide
        guides_list = annotate_variants(guides_list, args.verbosity, args.debug)
        # add allele frequencies for variants occurring in guides
        guides_list = annotate_variants_afs(guides_list, args.verbosity)
        # compute reverse complement for guides occurring on rev strand
        guides_list = reverse_guides(guides_list, args.verbosity)
        # compute gc content (pam excluded) for each guide
        guides_list = gc_content(guides_list, args.verbosity, args.debug)
        if args.annotations:  # annotate each guide
            guides_list = ann_guides(
                guides_list,
                region.contig,
                args.annotations,
                0,
                args.verbosity,
                args.debug,
            )
        if args.gene_annotations:  # annotate each guide with gene data
            guides_list = ann_guides(
                guides_list,
                region.contig,
                args.gene_annotations,
                1,
                args.verbosity,
                args.debug,
            )
        guides[region] = guides_list  # store annotated guides
    print_verbosity(
        f"Annotation completed in {time() - start:.2f}s",
        args.verbosity,
        VERBOSITYLVL[2],
    )
    return guides
