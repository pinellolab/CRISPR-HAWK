from crisprhawk import crisprme_data

import pytest


class DummyException(Exception):
    pass


def dummy_exception_handler(*args, **kwargs):
    raise DummyException("Exception handler called")


@pytest.fixture(autouse=True)
def patch_exception_handler(monkeypatch):
    monkeypatch.setattr(crisprme_data, "exception_handler", dummy_exception_handler)


def test_is_pam_right_upstream():
    header = "pam\tsgRNA_sequence"
    assert crisprme_data.is_pam_right(header) is True


def test_is_pam_right_downstream():
    header = "sgRNA_sequence\tpam"
    assert crisprme_data.is_pam_right(header) is False


def test_solve_pam_and_replacer():
    # IUPAC_ENCODER must contain 'R' for this test to work
    pamclass = "N[AG]G"
    # Simulate IUPAC_ENCODER for test
    orig_encoder = crisprme_data.IUPAC_ENCODER.copy()
    crisprme_data.IUPAC_ENCODER["AG"] = "R"
    result = crisprme_data.solve_pam(pamclass)
    assert result == "NRG"
    crisprme_data.IUPAC_ENCODER.clear()
    crisprme_data.IUPAC_ENCODER.update(orig_encoder)


def test_read_guides_report(tmp_path, monkeypatch):
    content = "pam\tsgRNA_sequence\tcol3\tcol4\tcol5\tcol6\nchr1\tA\tB\tC\tD\tE\nchr2\tF\tG\tH\tI\tJ\n"
    report_path = tmp_path / "report.txt"
    report_path.write_text(content)
    # Patch exception_handler to raise if called
    monkeypatch.setattr(crisprme_data, "exception_handler", dummy_exception_handler)
    right, fields = crisprme_data.read_guides_report(str(report_path), debug=False)
    assert right is True
    assert fields == [
        ["chr1", "A", "B", "C", "D", "E"],
        ["chr2", "F", "G", "H", "I", "J"],
    ]


def test_create_pam_file_and_create_guide_files(tmp_path, monkeypatch):
    outdir = tmp_path
    pamclass = "NGG"
    right = True
    guidelen = 20
    debug = False
    # Patch exception_handler to raise if called
    monkeypatch.setattr(crisprme_data, "exception_handler", dummy_exception_handler)
    crisprme_data.create_pam_file(pamclass, right, guidelen, str(outdir), debug)
    pamfile = outdir / "NGG.txt"
    assert pamfile.exists()
    content = pamfile.read_text().strip()
    assert content.startswith("NGG" + "N" * guidelen)
    # Test create_guide_files
    guides_data = [
        ["chr1", "A", "B", "NGG", "G" * guidelen, "NGG"],
        ["chr2", "C", "D", "NGG", "T" * guidelen, "NGG"],
    ]
    crisprme_data.create_guide_files(guides_data, right, str(outdir), debug)
    for g in guides_data:
        gfname = outdir / f"{g[0]}_{g[1]}_{g[2]}_{g[4]}_{g[3]}.txt"
        assert gfname.exists()
        gcontent = gfname.read_text().strip()
        assert gcontent.startswith("N" * len(g[3]) + g[4])
