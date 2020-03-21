import pytest

from aug.heredity.heredity import test_hwe


def test1():
    assert 0 == test_hwe(313, 102, 85)
