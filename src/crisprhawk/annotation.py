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
from .variant import adjust_multiallelic

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
    """Parses a variant identifier string into its chromosome, position, reference,
    and alternate alleles.

    This function splits the variant identifier, extracts the chromosome, position,
    reference, and alternate alleles, and normalizes them for multiallelic sites.

    Args:
        variant_id (str): The variant identifier string in the format 'chrom-pos-ref/alt'.

    Returns:
        Tuple[str, int, str, str]: A tuple containing the chromosome, normalized
            position, reference allele, and alternate allele.
    """
    parts = variant_id.split("-")
    chrom = parts[0]  # chromosome
    position = int(parts[1])  # variant position
    ref, alt = parts[2].split("/")
    # normalize for multiallelic sites
    ref_, alt_, position_ = adjust_multiallelic(ref, alt, position)
    return chrom, position_, ref_, alt_


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


def retrieve_guide_variants(guide: Guide) -> Set[str]:
    """Retrieves the set of variant identifiers associated with a guide.

    This function splits the guide's variants string into a set of individual
    variant identifiers.

    Args:
        guide (Guide): The guide object from which to retrieve variant identifiers.

    Returns:
        Set[str]: A set of variant identifier strings associated with the guide.
    """
    return set(guide.variants.split(","))  # split variants string


def is_reference_guide(guide_variants: Set[str], debug: bool) -> bool:
    """Determines if a guide is a reference guide based on its variant identifiers.

    This function checks if the guide's variants set contains only 'NA', indicating
    a reference guide. If 'NA' is present with other variants, it raises an exception.

    Args:
        guide_variants (Set[str]): Set of variant identifier strings for the guide.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        bool: True if the guide is a reference guide, False otherwise.

    Raises:
        ValueError: If 'NA' is present with other variants in the set.
    """
    if "NA" in guide_variants:  # reference guide
        if len(guide_variants) != 1:
            exception_handler(ValueError, "Forbidden NA variant", os.EX_DATAERR, debug)
        return True
    return False


def _create_variants_map(variants: Set[str]) -> List[Tuple[int, str]]:
    """Creates a mapping of variant positions to variant identifiers.

    This function parses each variant identifier to extract its genomic position
    and returns a list of (position, variant) tuples.

    Args:
        variants (Set[str]): A set of variant identifier strings.

    Returns:
        List[Tuple[int, str]]: A list of tuples mapping positions to variant
            identifiers.
    """
    variants_map = []  # variants map
    for variant in variants:
        _, pos, _, _ = _parse_variant(variant)
        variants_map.append((int(pos), variant))
    return variants_map


def search_variant(
    variants_map: List[Tuple[int, str]], pos: int, debug: bool
) -> Set[str]:
    """Searches for variant identifiers at a specific position in a variants map.

    This function returns a set of variant identifiers for the given position if found,
    otherwise raises an exception if the position is not present in the map.

    Args:
        variants_map (List[Tuple[int, str]]): List of (position, variant) tuples.
        pos (int): The genomic position to search for.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        Set[str]: A set of variant identifiers at the specified position.

    Raises:
        CrisprHawkAnnotationError: If the position is not found in the variants map.
    """
    variants = {v for p, v in variants_map if p == pos}
    if not variants:
        exception_handler(
            CrisprHawkAnnotationError,
            "Variant annotation requested for variant out of bounds",
            os.EX_DATAERR,
            debug,
        )
    return variants


def _find_insertion_stop(guideseq: str, debug: bool) -> int:
    """Finds the stop index of an insertion in a guide sequence.

    This function asserts that the insertion overlaps the start of the guide and
    returns the index of the first uppercase (reference) nucleotide, indicating
    where the insertion ends.

    Args:
        guideseq (str): The guide sequence containing the insertion.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        int: The index of the first reference nucleotide after the insertion.
    """
    # if insertion overlapping start, must not start with reference nts
    assert not guideseq[0].isupper() and not all(nt.isupper() for nt in guideseq)
    return next((i for i, nt in enumerate(guideseq) if nt.isupper()), 0)


def _check_insertion(
    guideseq: str, alt: str, posrel: int, pos: int, stop: int, is_snv: bool, debug: bool
) -> bool:
    """Checks if a guide sequence matches an insertion variant at a given position.

    This function determines whether the guide sequence context and the alternate allele
    are consistent with an insertion event, considering position and variant type.

    Args:
        guideseq (str): The guide sequence segment to check.
        alt (str): The alternate allele sequence.
        posrel (int): The relative position in the guide sequence.
        pos (int): The genomic position of the variant.
        stop (int): The stop position of the guide.
        is_snv (bool): Whether the variant is a single nucleotide variant.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        bool: True if the guide sequence matches the insertion variant, False otherwise.
    """
    if is_snv:  # do not even start
        return False
    if posrel == 0 and alt.endswith(
        guideseq.upper()[: _find_insertion_stop(guideseq, debug)]
    ):
        return True
    return bool(pos == stop and alt.startswith(guideseq.upper()))


def _check_snv(guideseq: str, alt: str) -> bool:
    """Checks if a guide sequence segment matches a single nucleotide variant (SNV)
    allele.

    This function returns True if the guide sequence is entirely lowercase and its
    uppercase form matches the alternate allele.

    Args:
        guideseq (str): The guide sequence segment to check.
        alt (str): The alternate allele sequence.

    Returns:
        bool: True if the guide sequence matches the SNV allele, False otherwise.
    """
    return guideseq.islower() and guideseq.upper() == alt


def polish_guide_variants(guide: Guide, variants: Set[str], debug: bool) -> str:
    """Validates and filters variants for a guide based on sequence and position.

    This function checks each variant for correct sequence context within the guide,
    handling indels and SNVs, and returns a comma-separated string of validated
    variant IDs.

    Args:
        guide (Guide): The guide object whose sequence and position map are used
            for validation.
        variants (Set[str]): Set of variant identifiers to validate.
        debug (bool): Flag to enable debug mode for error handling.

    Returns:
        str: A comma-separated string of validated variant identifiers for the guide.
    """
    variants_polished = set()  # polished set of variants for guide
    variants_map = _create_variants_map(variants)  # variants map
    variants_positions = {p for (p, _) in variants_map}
    variants = set(sorted(variants))
    for i, _ in enumerate(guide.guidepam):
        offset = 0  # offset for indels allele check
        if (p := guide.posmap[i]) in variants_positions:
            variant_ids = search_variant(variants_map, p, debug)
            for variant_id in variant_ids:
                # recover ref/alt alleles
                _, _, ref, alt = _parse_variant(variant_id)
                # indel, compute offset to retrieve allele
                if not (is_snv := _is_snv(ref, alt)):
                    offset = abs(len(ref) - len(alt)) if len(ref) < len(alt) else 0
                guide_segment = guide.guidepam[i : i + offset + 1]
                if _check_insertion(
                    guide_segment, alt, i, p, guide.stop, is_snv, debug
                ) or _check_snv(guide_segment, alt):
                    variants_polished.add(variant_id)
    return ",".join(sorted(variants_polished))


def annotate_variants(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    """Annotates guides with validated variant information based on their sequence
    context.

    This function processes each guide, determines if it is a reference or
    alternative guide, and updates its variants property with either 'NA' or a
    polished, validated set of variants.

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
        guide_variants = retrieve_guide_variants(guide)  # retrieve guide's variants
        if is_reference_guide(guide_variants, debug):  # reference guide
            guide.variants = "NA"
        else:  # alternative guide
            guide.variants = polish_guide_variants(guide, guide_variants, debug)
        guides_lst.append(guide)
    print_verbosity(
        f"Variants annotated in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides_lst

def _format_af(af: float) -> str:
    """Formats an allele frequency value for annotation output.

    This function returns the allele frequency as a string, using scientific 
    notation if it has more than four decimal digits.

    Args:
        af (float): The allele frequency value to format.

    Returns:
        str: The formatted allele frequency string.
    """
    assert isinstance(af, float)  # must be float
    s = f"{af:.10f}".rstrip("0").rstrip(".")  # check number of decimals
    decimal_digits = len(s.split(".")) if "." in s else 0
    return f"{af:.6e}" if decimal_digits > 3 else f"{af:.6f}"


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
                _format_af(guide.afs[v]) if str(guide.afs[v]) != "nan" else "NA"
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
