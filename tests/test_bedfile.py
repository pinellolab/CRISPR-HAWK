from crisprhawk.bedfile import (
    Bed,
    BedIterator,
    BedAnnotation,
    _parse_bed_line,
    _find_tbi,
)
from crisprhawk.coordinate import Coordinate

from tempfile import NamedTemporaryFile

import pytest
import os


class DummyFasta:
    def __init__(self, seq="ACTG"):
        self._seq = seq

    def fetch(self, coord):
        return self._seq


def make_bed_file(contents):
    f = NamedTemporaryFile(delete=False, mode="w", suffix=".bed")
    f.write(contents)
    f.close()
    return f.name


def test_bed_init_and_len(tmp_path):
    bed_content = "chr1\t10\t20\nchr2\t30\t40\n"
    bedfile = make_bed_file(bed_content)
    bed = Bed(bedfile, padding=0, debug=False)
    assert len(bed) == 2
    os.unlink(bedfile)


def test_bed_repr(tmp_path):
    bed_content = "chr1\t10\t20\n"
    bedfile = make_bed_file(bed_content)
    bed = Bed(bedfile, padding=0, debug=False)
    assert "stored regions=1" in repr(bed)
    os.unlink(bedfile)


def test_bed_getitem(tmp_path):
    bed_content = "chr1\t10\t20\nchr2\t30\t40\n"
    bedfile = make_bed_file(bed_content)
    bed = Bed(bedfile, padding=0, debug=False)
    coord = bed[0]
    assert isinstance(coord, Coordinate)
    coords = bed[:2]
    assert isinstance(coords, list)
    os.unlink(bedfile)


def test_bed_iter(tmp_path):
    bed_content = "chr1\t10\t20\nchr2\t30\t40\n"
    bedfile = make_bed_file(bed_content)
    bed = Bed(bedfile, padding=0, debug=False)
    items = list(bed)
    assert len(items) == 2
    os.unlink(bedfile)


def test_bed_extract_regions(tmp_path):
    bed_content = "chr1\t10\t20\n"
    bedfile = make_bed_file(bed_content)
    bed = Bed(bedfile, padding=0, debug=False)
    fastas = {"chr1": DummyFasta("ACTGACTGACTG")}
    regions = bed.extract_regions(fastas)
    assert hasattr(regions, "__iter__")
    os.unlink(bedfile)


def test_parse_bed_line_valid():
    line = "chr1\t10\t20"
    coord = _parse_bed_line(line, 1, 0, False)
    assert isinstance(coord, Coordinate)
    assert coord.contig == "chr1"
    assert coord.start == 10
    assert coord.stop == 20


def test_parse_bed_line_invalid_columns():
    with pytest.raises(ValueError):
        _parse_bed_line("chr1\t10", 1, 0, True)


def test_parse_bed_line_invalid_start_stop():
    with pytest.raises(TypeError):
        _parse_bed_line("chr1\ta\t20", 1, 0, True)


def test_parse_bed_line_stop_less_than_start():
    with pytest.raises(ValueError):
        _parse_bed_line("chr1\t20\t10", 1, 0, True)


def test_find_tbi(tmp_path):
    bedfile = tmp_path / "test.bed"
    bedfile.write_text("chr1\t10\t20\n")
    tbi_file = str(bedfile) + ".tbi"
    with open(tbi_file, "w") as f:
        f.write("dummy")
    assert _find_tbi(str(bedfile)) is True
    os.remove(tbi_file)
    assert _find_tbi(str(bedfile)) is False
