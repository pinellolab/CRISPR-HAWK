""" """

from .crisprhawk_error import (
    CrisprHawkCfdScoreError,
    CrisprHawkAzimuthScoreError,
    CrisprHawkRs3ScoreError,
    CrisprHawkDeepCpf1ScoreError,
    CrisprHawkAnnotationError,
)
from .exception_handlers import exception_handler
from .scores import azimuth, cfdon, rs3, deepcpf1
from .bedfile import BedAnnotation
from .guide import Guide, GUIDESEQPAD
from .utils import print_verbosity, flatten_list, VERBOSITYLVL
from .region import Region
from .pam import PAM, CASX, CPF1, SACAS9, SPCAS9, XCAS9
from .offtargets import search_offtargets

from collections import defaultdict
from typing import List, Dict, Union, Set
from time import time

import numpy as np

import os

ANNDIR = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "annotations"
)  # annotation data directory


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
        guides[i].set_azimuth_score(score)  # assign score to each guide
    print_verbosity(
        f"Azimuth scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def rs3_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
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
        guides[i].set_rs3_score(score)  # assign score to each guide
    print_verbosity(
        f"RS3 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def deepcpf1_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
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
        guides[i].set_deepcpf1_score(score)  # assign score to each guide
    print_verbosity(
        f"DeepCpf1 scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def group_guides_position(
    guides: List[Guide], debug: bool
) -> Dict[str, Dict[int, Union[Guide, List[Guide]]]]:
    # dictionary to map guides to positions (0 -> ref; 1 -> alt)
    pos_guide = defaultdict(lambda: {0: None, 1: []})
    for guide in guides:
        poskey = f"{guide.start}_{guide.strand}"
        if guide.samples == "REF":  # reference guide
            if pos_guide[poskey][0] is not None:
                exception_handler(
                    CrisprHawkCfdScoreError,
                    f"Duplicate REF guide at position {guide.start}? CFDon calculation failed",
                    os.EX_DATAERR,
                    debug,
                )
            pos_guide[poskey][0] = guide  # type: ignore
        pos_guide[poskey][1].append(guide)  # type: ignore
    return pos_guide  # type: ignore


def cfdon_score(guides: List[Guide], verbosity: int, debug: bool) -> List[Guide]:
    print_verbosity("Computing CFDon score", verbosity, VERBOSITYLVL[3])
    start = time()  # cfdon start time
    guide_groups = group_guides_position(guides, debug)  # group guides by positions
    for _, gg in guide_groups.items():
        try:
            cfdon_scores = cfdon(gg[0], gg[1], debug)  # type: ignore
        except Exception as e:
            exception_handler(
                CrisprHawkCfdScoreError,
                "CFDon score calculation failed",
                os.EX_DATAERR,
                debug,
                e,
            )
        for i, score in enumerate(cfdon_scores):
            gg[1][i].set_cfdon_score(score)  # type: ignore
    # revert grouped guides by position into list
    guides = flatten_list([gg[1] for _, gg in guide_groups.items()])  # type: ignore
    print_verbosity(
        f"CFDon scores computed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides


def polish_variants_annotation(guide: Guide, variants: Set[str]) -> Set[str]:
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
            guide.set_variants(",".join(sorted(guide_vars)))
            guides_lst.append(guide)
    print_verbosity(
        f"Variants annotated in {time() - start:.2f}s", verbosity, VERBOSITYLVL[3]
    )
    return guides_lst


def annotate_variants_afs(
    guides: List[Guide], verbosity: int, debug: bool
) -> List[Guide]:
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
        guide.set_allele_freqs(afs)
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
    guide.set_func_ann(annotation)
    return guide


def _retrieve_gene_name(field: str) -> str:
    i = field.find("gene_name=")
    if i == -1:  # gene name not found
        return ""
    j = field.find(";", i + 10)
    return field[i + 10 : j]  # return gene name


def _geneann(
    guide: Guide, bedannotation: BedAnnotation, contig: str, debug: bool
) -> Guide:
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
    guide.set_gene_ann(annotation)
    return guide


def ann_guides(
    guides: List[Guide],
    contig: str,
    annotations: List[str],
    atype: int,  # 0 -> regular annotation; 1 -> gene annotation
    verbosity: int,
    debug: bool,
) -> List[Guide]:
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
    return guides_ann if guides_ann else guides


def annotate_guides(
    guides: Dict[Region, List[Guide]],
    annotations: List[str],
    gene_annotations: List[str],
    pam: PAM,
    genome: str,
    estimate_offtargets: bool,
    outdir: str,
    threads: int,
    verbosity: int,
    debug: bool,
) -> Dict[Region, List[Guide]]:
    # annotate guides with scores, variants and adjust positions
    print_verbosity("Annotating guides", verbosity, VERBOSITYLVL[1])
    start = time()  # annotation start time
    for region, guides_list in guides.items():
        # set variants for current guide
        guides_list = annotate_variants(guides_list, verbosity, debug)
        # add allele frequencies for variants occurring in guides
        guides_list = annotate_variants_afs(guides_list, verbosity, debug)
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
        # TODO: use crispritz instead
        if estimate_offtargets:  # estimate off-targets for each guide
            guides_list = search_offtargets(
                guides_list,
                pam,
                genome,
                region,
                "",  # TODO: fix
                "",
                threads,
                outdir,
                verbosity,
                debug,
            )
        # guides[region] = guides_list  # store annotated guides
    print_verbosity(
        f"Annotation completed in {time() - start:.2f}s", verbosity, VERBOSITYLVL[2]
    )
    return guides
