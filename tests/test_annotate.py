from unittest.mock import MagicMock, patch
from crisprhawk.annotate import reverse_guides, azimuth_score, rs3_score, deepcpf1_score, group_guides_position
from crisprhawk.utils import prepare_package

import numpy as np

import pytest


@pytest.fixture(scope="session")
def uncompress_data():
    prepare_package()


class DummyGuide:
    def __init__(self, seq, strand, samples):
        self.sequence = seq
        self.strand = strand
        self.samples = samples
        self.reverse_complement_called = False

    def reverse_complement(self):
        self.reverse_complement_called = True
        # For test, just reverse the string and swap case to simulate change
        self.sequence = self.sequence[::-1].swapcase()

    def set_azimuth_score(self, score):
        self.azimuth_score = score

    def set_rs3_score(self, score):
        self.rs3_score = score

    def set_deepcpf1_score(self, score):
        self.deepcpf1_score = score

@pytest.mark.parametrize(
    "guides,verbosity,expected_sequences,expected_rc_calls",
    [
        # No guides
        ([], 1, [], [],),
        # All guides on forward strand
        (
            [DummyGuide("ATGC", 0, ""), DummyGuide("GGCC", 0, "")],
            2,
            ["ATGC", "GGCC"],
            [False, False],
        ),
        # All guides on reverse strand
        (
            [DummyGuide("ATGC", 1, ""), DummyGuide("GGCC", 1, "")],
            3,
            ["CGTA", "CCGG"],  # reversed and swapcase
            [True, True],
        ),
        # Mixed strands
        (
            [DummyGuide("ATGC", 0, ""), DummyGuide("GGCC", 1, "")],
            0,
            ["ATGC", "CCGG"],
            [False, True],
        ),
        # Single guide, reverse strand
        (
            [DummyGuide("TTAA", 1, "")],
            1,
            ["AATT"],
            [True],
        ),
        # Single guide, forward strand
        (
            [DummyGuide("TTAA", 0, "")],
            1,
            ["TTAA"],
            [False],
        ),
    ],
    ids=[
        "empty-list",
        "all-forward",
        "all-reverse",
        "mixed-strands",
        "single-reverse",
        "single-forward",
    ]
)
def test_reverse_guides_happy_and_edge_cases(guides, verbosity, expected_sequences, expected_rc_calls):
    # Arrange
    with patch("crisprhawk.annotate.print_verbosity") as mock_print, \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]):
        # Act
        result = reverse_guides(guides, verbosity)

        # Assert
        assert result == guides
        for guide, expected_seq, expected_rc in zip(result, expected_sequences, expected_rc_calls):
            assert guide.sequence.upper() == expected_seq
            assert guide.reverse_complement_called == expected_rc  # pyright: ignore[reportAttributeAccessIssue]
        # print_verbosity should be called twice (start and end)
        assert mock_print.call_count == 2

@pytest.mark.parametrize(
    "guides,verbosity",
    [
        ([DummyGuide("ATGC", 1, "")], 1),
        ([DummyGuide("ATGC", 0, "")], 1),
    ],
    ids=["reverse-guide-print-error", "forward-guide-print-error"]
)
def test_reverse_guides_print_verbosity_error(monkeypatch, guides, verbosity):
    # Arrange
    def raise_on_print(*args, **kwargs):
        raise RuntimeError("print_verbosity failed")
    monkeypatch.setattr("crisprhawk.annotate.print_verbosity", raise_on_print)
    monkeypatch.setattr("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3])

    # Act & Assert
    with pytest.raises(RuntimeError, match="print_verbosity failed"):
        reverse_guides(guides, verbosity)

def test_reverse_guides_reverse_complement_error(monkeypatch):
    # Arrange
    class BadGuide(DummyGuide):
        def reverse_complement(self):
            raise ValueError("reverse_complement failed")
    guides = [BadGuide("ATGC", 1, "")]
    monkeypatch.setattr("crisprhawk.annotate.print_verbosity", lambda *a, **k: None)
    monkeypatch.setattr("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3])

    # Act & Assert
    with pytest.raises(ValueError, match="reverse_complement failed"):
        reverse_guides(guides, 1)  # type: ignore


@pytest.mark.parametrize(
    "guides,verbosity,debug,expected_scores,expected_return,case_id",
    [
        # No guides
        ([], 1, False, [], [], "empty-list"),
        # Single guide, happy path
        ([DummyGuide("AAAACCCCGGGGTTTT", 0, "")], 2, False, [0.5], [0.5], "single-guide"),
        # Multiple guides, happy path
        (
            [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 0, "")],
            3,
            False,
            [0.1, 0.2],
            [0.1, 0.2],
            "multiple-guides"
        ),
    ],
    ids=lambda x: x if isinstance(x, str) else None
)
def test_azimuth_score_happy_and_edge_cases(guides, verbosity, debug, expected_scores, expected_return, case_id):
    # Arrange
    with patch("crisprhawk.annotate.print_verbosity") as mock_print, \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.azimuth") as mock_azimuth:
        if guides:
            mock_azimuth.return_value = np.array(expected_scores)
        # Act
        result = azimuth_score(guides, verbosity, debug)
        # Assert
        assert result == guides
        for guide, expected in zip(result, expected_return):
            assert guide.azimuth_score == expected
        if guides:
            assert mock_print.call_count == 2
        else:
            assert mock_print.call_count == 0

def test_azimuth_score_assertion_error():
    # Arrange
    guides = [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 0, "")]
    with patch("crisprhawk.annotate.print_verbosity"), \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.azimuth", return_value=np.array([0.1])):  # wrong length
        # Act & Assert
        with pytest.raises(AssertionError):
            azimuth_score(guides, 1, False) # pyright: ignore[reportArgumentType]


    

@pytest.mark.parametrize(
    "guides,verbosity,debug,rs3_return,expected_scores,expected_return,case_id",
    [
        # No guides
        ([], 1, False, [], [], [], "empty-list"),
        # Single guide, happy path
        ([DummyGuide("AAAACCCCGGGGTTTT", 0, "")], 2, False, [0.5], [0.5], [0.5], "single-guide"),
        # Multiple guides, happy path
        (
            [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 0, "")],
            3,
            False,
            [0.1, 0.2],
            [0.1, 0.2],
            [0.1, 0.2],
            "multiple-guides"
        ),
    ],
    ids=lambda x: x if isinstance(x, str) else None
)
def test_rs3_score_happy_and_edge_cases(guides, verbosity, debug, rs3_return, expected_scores, expected_return, case_id):
    # Arrange
    with patch("crisprhawk.annotate.print_verbosity") as mock_print, \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.rs3") as mock_rs3:
        if guides:
            mock_rs3.return_value = expected_scores
        # Act
        result = rs3_score(guides, verbosity, debug)
        # Assert
        assert result == guides
        for guide, expected in zip(result, expected_return):
            assert guide.rs3_score == expected
        if guides:
            assert mock_print.call_count == 2
        else:
            assert mock_print.call_count == 0

def test_rs3_score_assertion_error():
    # Arrange
    guides = [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 1, "")]
    with patch("crisprhawk.annotate.print_verbosity"), \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.rs3", return_value=[0.1]):  # wrong length
        # Act & Assert
        with pytest.raises(AssertionError):
            rs3_score(guides, 1, False) # pyright: ignore[reportArgumentType]

@pytest.mark.parametrize(
    "guides,verbosity,debug,deepcpf1_return,expected_scores,expected_return,case_id",
    [
        # No guides
        ([], 1, False, [], [], [], "empty-list"),
        # Single guide, happy path
        ([DummyGuide("AAAACCCCGGGGTTTT", 1, "")], 2, False, [0.5], [0.5], [0.5], "single-guide"),
        # Multiple guides, happy path
        (
            [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 1, "")],
            3,
            False,
            [0.1, 0.2],
            [0.1, 0.2],
            [0.1, 0.2],
            "multiple-guides"
        ),
    ],
    ids=lambda x: x if isinstance(x, str) else None
)
def test_deepcpf1_score_happy_and_edge_cases(guides, verbosity, debug, deepcpf1_return, expected_scores, expected_return, case_id):
    # Arrange
    with patch("crisprhawk.annotate.print_verbosity") as mock_print, \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.deepcpf1") as mock_deepcpf1:
        if guides:
            mock_deepcpf1.return_value = expected_scores
        # Act
        result = deepcpf1_score(guides, verbosity, debug)
        # Assert
        assert result == guides
        for guide, expected in zip(result, expected_return):
            assert guide.deepcpf1_score == expected
        if guides:
            assert mock_print.call_count == 2
        else:
            assert mock_print.call_count == 0


def test_deepcpf1_score_assertion_error():
    # Arrange
    guides = [DummyGuide("AAAACCCCGGGGTTTT", 0, ""), DummyGuide("TTTTGGGGCCCCAAAA", 0, "")]
    with patch("crisprhawk.annotate.print_verbosity"), \
         patch("crisprhawk.annotate.VERBOSITYLVL", [0, 1, 2, 3]), \
         patch("crisprhawk.annotate.GUIDESEQPAD", 4), \
         patch("crisprhawk.annotate.deepcpf1", return_value=[0.1]):  # wrong length
        # Act & Assert
        with pytest.raises(AssertionError):
            deepcpf1_score(guides, 1, False) # pyright: ignore[reportArgumentType]


class DummyGuide2:
    def __init__(self, start, strand, samples):
        self.start = start
        self.strand = strand
        self.samples = samples


@pytest.mark.parametrize(
    "guides,debug,expected_keys,expected_refs,expected_guides,case_id",
    [
        # No guides
        ([], False, [], [], [], "empty-list"),
        # Single guide, not reference
        ([DummyGuide2(100, 0, "ALT")], False, ["100_0"], [None], [[DummyGuide2(100, 0, "ALT")]], "single-nonref"),
        # Single guide, reference
        ([DummyGuide2(200, 1, "REF")], False, ["200_1"], [DummyGuide2(200, 1, "REF")], [[DummyGuide2(200, 1, "REF")]], "single-ref"),
        # Multiple guides, same position, one reference, one alt
        ([DummyGuide2(300, 0, "REF"), DummyGuide2(300, 0, "ALT")], False, ["300_0"], [DummyGuide2(300, 0, "REF")], [[DummyGuide2(300, 0, "REF"), DummyGuide2(300, 0, "ALT")]], "multi-same-pos-ref-alt"),
        # Multiple guides, different positions and strands
        (
            [DummyGuide2(400, 0, "REF"), DummyGuide2(401, 1, "ALT"), DummyGuide2(400, 0, "ALT")],
            False,
            ["400_0", "401_1"],
            [DummyGuide2(400, 0, "REF"), None],
            [[DummyGuide2(400, 0, "REF"), DummyGuide2(400, 0, "ALT")], [DummyGuide2(401, 1, "ALT")]],
            "multi-diff-pos"
        ),
    ],
    ids=lambda x: x if isinstance(x, str) else None
)
def test_group_guides_position_happy_and_edge_cases(guides, debug, expected_keys, expected_refs, expected_guides, case_id):
    # Arrange
    # Act
    result = group_guides_position(guides, debug)
    # Assert
    assert set(result.keys()) == set(expected_keys)
    for key, expected_ref, expected_guides_list in zip(expected_keys, expected_refs, expected_guides):
        ref, guides_list = result[key]
        if expected_ref is None:
            assert ref is None
        else:
            assert ref.start == expected_ref.start and ref.strand == expected_ref.strand and ref.samples == expected_ref.samples # type: ignore
        assert len(guides_list) == len(expected_guides_list)
        for g, eg in zip(guides_list, expected_guides_list):
            assert g.start == eg.start and g.strand == eg.strand and g.samples == eg.samples

def test_group_guides_position_duplicate_ref_triggers_exception_handler():
    # Arrange
    guides = [DummyGuide2(500, 0, "REF"), DummyGuide2(500, 0, "REF")]
    with patch("crisprhawk.annotate.exception_handler") as mock_handler, \
         patch("crisprhawk.annotate.CrisprHawkCfdScoreError", new=Exception):
        # Act
        group_guides_position(guides, debug=True) # type: ignore
        # Assert
        assert mock_handler.call_count == 1
        args, kwargs = mock_handler.call_args
        assert args[0] == Exception
        assert "Duplicate REF guide at position 500" in args[1]
        assert args[2] == 65  # os.EX_DATAERR
        assert args[3] is True

