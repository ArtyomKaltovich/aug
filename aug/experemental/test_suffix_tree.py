import pytest
import random

from aug.experemental.SuffixTree import SuffixTree, _TreeElem, SuffixTreeBuildingMethod


@pytest.fixture(scope='module')
def seed():
    seed = random.randrange(0, 10 ** 5)
    print(f"seed={seed}")
    return random.seed(seed)


def test_extract_general_prefix():
    tree = _TreeElem("a")
    assert 1 == tree.extract_general_prefix(_TreeElem("anna"))
    tree = _TreeElem("na")
    assert 1 == tree.extract_general_prefix(_TreeElem("nna"))


def test_one_string1():
    expected = "Elem(node=None, children=" \
               "Elem(node=a, children=Elem(node=na, children=Elem(node=na, children=Elem(node=1, children=)), Elem(node=1, children=)), Elem(node=1, children=)), " \
               "Elem(node=banana, children=Elem(node=1, children=)), " \
               "Elem(node=na, children=Elem(node=na, children=Elem(node=1, children=)), Elem(node=1, children=)))"
    assert expected == SuffixTree("banana", method=SuffixTreeBuildingMethod.naive).to_graphviz()


def test_one_string2():
    expected = "Elem(node=None, children=" \
               "Elem(node=a, children=Elem(node=nna, children=Elem(node=1, children=)), " \
               "Elem(node=1, children=)), " \
               "Elem(node=n, children=Elem(node=a, children=Elem(node=1, children=)), Elem(node=na, children=Elem(node=1, children=))))"
    assert expected == str(SuffixTree("anna", method=SuffixTreeBuildingMethod.naive))


def test_one_string3():
    expected = "Elem(node=None, children=" \
               "Elem(node=a, children=Elem(node=na, children=Elem(node=na, children=Elem(node=1, children=)), Elem(node=1, children=)), Elem(node=1, children=)), " \
               "Elem(node=banana, children=Elem(node=1, children=)), " \
               "Elem(node=na, children=Elem(node=na, children=Elem(node=1, children=)), Elem(node=1, children=)))"
    assert expected == str(SuffixTree("banana", method=SuffixTreeBuildingMethod.naive))


def test_one_letter_string():
    expected = "Elem(node=None, " \
               "children=Elem(node=a, " \
                    "children=Elem(node=a, " \
                        "children=Elem(node=a, " \
                            "children=Elem(node=a, children=Elem(node=1, children=)), " \
                        "Elem(node=1, children=)), " \
                    "Elem(node=1, children=)), " \
               "Elem(node=1, children=)))"
    assert expected == str(SuffixTree("aaaa", method=SuffixTreeBuildingMethod.naive))


def test_two_string():
    expected = "Elem(node=None, " \
               "children=Elem(node=a, children=Elem(node=na, children=Elem(node=na, children=Elem(node=s, children=Elem(node=2, children=)), " \
                    "Elem(node=1, children=)), Elem(node=s, children=Elem(node=2, children=)), Elem(node=1, children=)), Elem(node=s, children=Elem(node=2, children=)), Elem(node=1, children=)), " \
               "Elem(node=banana, children=Elem(node=1, children=)), " \
               "Elem(node=na, children=Elem(node=na, children=Elem(node=s, children=Elem(node=2, children=)), Elem(node=1, children=)), Elem(node=s, children=Elem(node=2, children=)), Elem(node=1, children=)), Elem(node=s, children=Elem(node=2, children=)))"
    assert expected == str(SuffixTree(["banana", "ananas"], method=SuffixTreeBuildingMethod.naive))


def test_two_string2():
    expected = "Elem(node=None, children=" \
               "Elem(node=a, children=Elem(node=b, children=Elem(node=a, children=Elem(node=2, children=)), Elem(node=ba, children=Elem(node=1, children=))), Elem(node=1, children=), Elem(node=2, children=)), " \
               "Elem(node=b, children=Elem(node=a, children=Elem(node=ba, children=Elem(node=2, children=)), Elem(node=1, children=), Elem(node=2, children=)), Elem(node=ba, children=Elem(node=1, children=))))"
    assert expected == str(SuffixTree(["abba", "baba"], method=SuffixTreeBuildingMethod.naive))


def test_three_string():
    expected = "Elem(node=None, children=" \
               "Elem(node=a, children=Elem(node=b, children=Elem(node=a, children=Elem(node=b, children=Elem(node=3, children=)), " \
                    "Elem(node=2, children=)), Elem(node=ba, children=Elem(node=1, children=)), Elem(node=3, children=)), Elem(node=1, children=), Elem(node=2, children=)), " \
               "Elem(node=b, children=Elem(node=a, children=Elem(node=b, children=Elem(node=a, children=Elem(node=2, children=)), " \
                    "Elem(node=3, children=)), Elem(node=1, children=), Elem(node=2, children=)), Elem(node=ba, children=Elem(node=1, children=)), Elem(node=3, children=)))"
    assert expected == str(SuffixTree(["abba", "baba", "abab"], method=SuffixTreeBuildingMethod.naive))
