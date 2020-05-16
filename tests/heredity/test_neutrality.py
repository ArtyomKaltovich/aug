import pytest

from aug.heredity.heredity import selection_balance


class TestDNDS:
    @pytest.mark.parametrize("seq1", ["ATTAGTCGTTCATCC"])
    @pytest.mark.parametrize("seq2", ["ATTAGACGTTCCTCC"])
    @pytest.mark.parametrize("expected", [pytest.approx(0.318)])
    def test_simple(self, seq1, seq2, expected):
        assert expected == selection_balance(seq1, seq2)
