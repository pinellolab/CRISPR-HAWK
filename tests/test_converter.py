import os
import tempfile
import shutil
import pytest

from crisprhawk.converter import (
    tabix_index,
    load_vcf,
    format_ac,
    variant_observed,
    _update_header,
    _asses_genotype,
    _format_vrecord,
    _convert,
    _compress,
    convert_vcf,
    convert_gnomad_vcf,
    GNOMADPOPS,
    GT,
    GTLINE,
)
from crisprhawk.crisprhawk_error import CrisprHawkConverterError

import pysam


class DummyVariant:
    def __init__(
        self,
        info,
        chrom="1",
        pos=1,
        id=".",
        ref="A",
        alts=("C",),
        qual=100,
        filter_keys=("PASS",),
    ):
        self.info = info
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alts = alts
        self.qual = qual
        self.filter = {k: None for k in filter_keys}


@pytest.fixture(autouse=True)
def patch_exception_handler(monkeypatch):
    monkeypatch.setattr(
        "crisprhawk.converter.exception_handler",
        lambda exc_type, msg, code, debug, *args: (_ for _ in ()).throw(exc_type(msg)),
    )


def test_format_ac():
    acs = format_ac(joint=False)
    assert acs == [f"AC_{p}" for p in GNOMADPOPS]
    acs_joint = format_ac(joint=True)
    assert acs_joint == [f"AC_joint_{p}" for p in GNOMADPOPS]


def test_variant_observed():
    assert variant_observed((0, 0, 0)) is False
    assert variant_observed((0, 1, 0)) is True
    assert variant_observed((1,)) is True


def test_update_header_adds_gt_and_samples():
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    new_header = _update_header(header, joint=False)
    assert "GT" in new_header
    for pop in GNOMADPOPS:
        assert pop in new_header


def test_update_header_joint_replaces_af():
    header = pysam.VariantHeader()
    header.add_line(
        '##INFO=<ID=AF_joint,Number=A,Type=Float,Description="Allele Frequency">'
    )
    new_header = _update_header(header, joint=True)
    assert "<ID=AF," in new_header


def test_asses_genotype_success():
    variant = DummyVariant(info={f"AC_{p}": (0,) for p in GNOMADPOPS})
    ac_formatted = [f"AC_{p}" for p in GNOMADPOPS]
    gt = _asses_genotype(variant, ac_formatted, debug=False)
    assert gt == "\t".join([GT[0]] * len(GNOMADPOPS))
    variant2 = DummyVariant(info={f"AC_{p}": (1,) for p in GNOMADPOPS})
    gt2 = _asses_genotype(variant2, ac_formatted, debug=False)
    assert gt2 == "\t".join([GT[1]] * len(GNOMADPOPS))


def test_asses_genotype_error():
    variant = DummyVariant(info={})
    ac_formatted = [f"AC_{p}" for p in GNOMADPOPS]
    with pytest.raises(CrisprHawkConverterError):
        _asses_genotype(variant, ac_formatted, debug=True)


def test_format_vrecord_with_af():
    variant = DummyVariant(info={"AF": [0.1, 0.2]}, alts=("C", "G"))
    genotypes = "0/0\t0/1"
    rec = _format_vrecord(variant, genotypes)
    assert "AF=0.1,0.2" in rec
    assert rec.endswith(genotypes)


def test_format_vrecord_missing_af():
    variant = DummyVariant(info={}, alts=("C", "G"))
    genotypes = "0/0\t0/1"
    rec = _format_vrecord(variant, genotypes)
    assert "AF=0.0,0.0" in rec


def test_convert_and_compress(tmp_path, monkeypatch):
    # Create a dummy VCF file
    vcf_path = tmp_path / "test.vcf"
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1\t.\tA\tC\t100\tPASS\tAF=0.1\n")
    # Patch pysam.tabix_index and find_tbi
    monkeypatch.setattr("crisprhawk.converter.pysam.tabix_index", lambda *a, **k: None)
    monkeypatch.setattr("crisprhawk.converter.find_tbi", lambda fname: True)
    # tabix_index should return the tbi path
    idx = tabix_index(str(vcf_path), verbosity=0, debug=False)
    assert idx.endswith(".tbi")
    # load_vcf should return a VariantFile
    vfile = load_vcf(str(vcf_path), verbosity=0, debug=False)
    assert isinstance(vfile, pysam.VariantFile)


def test_convert_gnomad_vcf_parallel(monkeypatch, tmp_path):
    # Patch convert_vcf to just record calls
    called = []
    monkeypatch.setattr(
        "crisprhawk.converter.convert_vcf",
        lambda vcf_fname, **kwargs: called.append(vcf_fname),
    )
    monkeypatch.setattr(
        "crisprhawk.converter.multiprocessing.Pool",
        lambda processes: type(
            "Pool",
            (),
            {
                "map": lambda self, func, files: [func(f) for f in files],
                "close": lambda self: None,
                "join": lambda self: None,
            },
        )(),
    )
    files = [str(tmp_path / f"f{i}.vcf") for i in range(2)]
    convert_gnomad_vcf(
        files,
        joint=False,
        keep=True,
        suffix="s",
        outdir=str(tmp_path),
        threads=1,
        verbosity=0,
        debug=False,
    )
    assert called == files
