from crisprhawk import utils

import tempfile
import pytest
import shutil
import os


def test_reverse_complement_basic():
    assert utils.reverse_complement("ATGC", False) == "GCAT"
    assert utils.reverse_complement("atgc", False) == "gcat"
    assert utils.reverse_complement("N", False) == "N"

def test_reverse_complement_invalid():
    with pytest.raises(SystemExit):
        utils.reverse_complement("ATXG", False)

def test_warning_and_print_verbosity(capsys):
    utils.warning("test warning", 1)
    captured = capsys.readouterr()
    assert "WARNING: test warning." in captured.err

    utils.print_verbosity("test message", 2, 1)
    captured = capsys.readouterr()
    assert "test message" in captured.out

def test_adjust_guide_position():
    assert utils.adjust_guide_position(10, 5, 3, True) == 10
    assert utils.adjust_guide_position(10, 5, 3, False) == 5

def test_round_score():
    assert utils.round_score(1.234567) == 1.2346
    assert utils.round_score(1.0) == 1.0

def test_flatten_list():
    assert utils.flatten_list([[1, 2], [3], [], [4, 5]]) == [1, 2, 3, 4, 5]
    assert utils.flatten_list([]) == []

def test_match_iupac():
    assert utils.match_iupac("A", "A")
    assert utils.match_iupac("A", "N")
    assert not utils.match_iupac("A", "C")
    assert utils.match_iupac("AG", "RR")
    assert not utils.match_iupac("AG", "YY")
    assert not utils.match_iupac("AA", "AG")

def test_dna2rna():
    assert utils.dna2rna("ATGC") == "AUGC"
    assert utils.dna2rna("atgc") == "augc"
    assert utils.dna2rna("TTTT") == "UUUU"

def test_create_and_remove_folder():
    temp_dir = tempfile.mkdtemp()
    new_dir = os.path.join(temp_dir, "testdir")
    assert utils.create_folder(new_dir) == new_dir
    assert os.path.isdir(new_dir)
    utils.remove_folder(new_dir)
    assert not os.path.exists(new_dir)
    shutil.rmtree(temp_dir, ignore_errors=True)

def test_remove_file():
    temp_dir = tempfile.mkdtemp()
    file_path = os.path.join(temp_dir, "testfile.txt")
    with open(file_path, "w") as f:
        f.write("test")
    assert os.path.isfile(file_path)
    utils.remove_file(file_path)
    assert not os.path.exists(file_path)
    shutil.rmtree(temp_dir, ignore_errors=True)

def test_command_exists():
    assert utils.command_exists("python")
    assert not utils.command_exists("some_nonexistent_command_12345")

def test_is_lowercase():
    assert utils.is_lowercase("abc")
    assert not utils.is_lowercase("ABC")
    assert utils.is_lowercase("aBC")
    assert not utils.is_lowercase("")