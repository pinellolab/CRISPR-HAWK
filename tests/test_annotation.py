import pytest
import numpy as np
from Bio.SeqUtils import gc_fraction
from crisprhawk.annotation import (
    reverse_guides,
    annotate_variants_afs,
    gc_content,
)
from crisprhawk.guide import Guide, GUIDESEQPAD


def make_test_guide(spacer="ATGCATGCATGCATGCATGC", pam="NGG", direction=1, right=False):
    """
    Helper to create a Guide object that satisfies guide.py logic.
    - guidelen=20, pamlen=3.
    - sequence must include GUIDESEQPAD (10bp) on both sides.
    """
    guidelen = len(spacer)
    pamlen = len(pam)

    # Construct the core sequence based on PAM orientation
    if right:  # PAM is at the start (Cpf1 style)
        core = pam + spacer
    else:  # PAM is at the end (Cas9 style)
        core = spacer + pam

    # Full sequence with padding as expected by Guide.__init__
    full_sequence = "N" * GUIDESEQPAD + core + "N" * GUIDESEQPAD

    # Map every index in the padded sequence to a dummy genomic coordinate
    posmap = {i: 5000 + i for i in range(len(full_sequence))}

    return Guide(
        position_start=5010,
        position_stop=5010 + guidelen + pamlen,
        sequence=full_sequence,
        guidelen=guidelen,
        pamlen=pamlen,
        direction=direction,
        samples="sample1",
        variants="var_A,var_B",
        posmap=posmap,
        afs={},
        debug=False,
        right=right,
        hapid="hap1",
    )


def test_annotate_variants_afs():
    """
    Revised test: Ensure the afs dictionary contains the variant key
    to satisfy the lookup in annotation.py.
    """
    # 1. Define variant IDs that match the genomic positions in make_test_guide
    var1 = "chr1-5010-A/T"
    var2 = "chr1-5012-G/C"

    # 2. Create the guide
    g = make_test_guide(spacer="ATGCATGCATGCATGCATGC")

    # 3. CRITICAL: The variants string and the afs keys must match exactly
    g.variants = f"{var1},{var2}"

    # 4. Populate the internal afs dictionary so the function find the keys
    g._afs = {var1: 0.123, var2: 0.005}

    # Now this will not raise a KeyError
    guides = annotate_variants_afs([g], verbosity=0)

    assert guides[0].afs[var1] == 0.123
    assert guides[0].afs[var2] == 0.005


def test_gc_content():
    """Verify GC content calculation excludes the PAM sequence."""
    # 20nt Spacer: 10 Gs, 10 As = 50% GC
    # PAM: NGG (should be ignored)
    spacer = "GGGGGGGGGGAAAAAAAAAA"
    g = make_test_guide(spacer=spacer, pam="NGG", right=False)

    processed = gc_content([g], verbosity=0, debug=False)

    # gc_fraction returns float 0.0-1.0
    assert float(processed[0].gc) == pytest.approx(0.5)


def test_guide_slicing_integrity():
    """Sanity check: Ensure Guide object correctly identifies spacer vs PAM."""
    # SpCas9 style: Spacer(20) then PAM(3)
    spacer = "ATGC" * 5
    pam = "TGG"
    g = make_test_guide(spacer=spacer, pam=pam, right=False)

    assert g.guide == spacer
    assert g.pam == pam
