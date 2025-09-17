from crisprhawk.coordinate import Coordinate

import pytest


def test_coordinate_init_and_properties():
    c = Coordinate("chr1", 100, 200, 10)
    assert c.contig == "chr1"
    assert c.start == 90  # start - padding
    assert c.stop == 210  # stop + padding
    assert c.startp == 100
    assert c.stopp == 200
    assert c._padding == 10


def test_coordinate_init_padding_does_not_go_below_zero():
    c = Coordinate("chr2", 5, 10, 10)
    assert c.start == 0  # start - padding, but not below 0
    assert c.stop == 20


def test_coordinate_init_stop_less_than_start():
    with pytest.raises(ValueError):
        Coordinate("chr1", 100, 50, 5)


def test_coordinate_str_and_repr():
    c = Coordinate("chr3", 100, 200, 5)
    assert str(c) == "chr3:100-200"
    # repr uses start+padding and stop-padding
    # but if padding=5, repr will show 105-195
    assert "chr3:105-195" in repr(c)
    assert "padding=5" in repr(c)


def test_coordinate_equality():
    c1 = Coordinate("chr1", 100, 200, 10)
    c2 = Coordinate("chr1", 100, 200, 10)
    c3 = Coordinate("chr1", 100, 201, 10)
    c4 = Coordinate("chr2", 100, 200, 10)
    assert c1 == c2
    assert not (c1 == c3)
    assert not (c1 == c4)
    assert c1 != "not_a_coordinate"
