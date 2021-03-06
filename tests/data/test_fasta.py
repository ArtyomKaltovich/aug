import pytest

from aug.seq.seq import gc_rate, profile
from aug.data.fasta import *
from tests import FLOAT_EQUALITY_ACCURACY


def test_fasta_file_do(base_data_path):
    def helper(id, string):
        return id, string

    actual = fasta_data_do(base_data_path + r"test_fasta.txt", helper, [])
    expected = [("Rosalind_6404", "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("Rosalind_5959", "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("Rosalind_0808", "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == actual


def test_read_fasta(base_data_path):
    expected = [("Rosalind_6404", "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("Rosalind_5959", "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("Rosalind_0808", "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == read_fasta(base_data_path + r"test_fasta.txt")


def test_read_fasta_without_id(base_data_path):
    expected = [("CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == read_fasta(base_data_path + r"test_fasta.txt", without_id=True)


def test_fasta_file_do_max(base_data_path):
    id_max = None
    gc_max = None

    def helper(id, string):
        nonlocal id_max, gc_max
        gc = gc_rate(string, percent=True)
        if not gc_max or gc > gc_max:
            id_max = id
            gc_max = gc

    fasta_data_do(base_data_path + r"test_fasta.txt", helper)
    assert "Rosalind_0808" == id_max
    assert 60.919540 == pytest.approx(gc_max)


def test_fasta_profile(base_data_path):
    profile_matrix = None
    for id, dna in fasta_file_iter(base_data_path + "test_profile_fasta.txt"):
        profile_matrix = profile(dna, profile_matrix)
    expected = {'A': [5, 1, 0, 0, 5, 5, 0, 0], 'C': [0, 0, 1, 4, 2, 0, 6, 1], 'G': [1, 1, 6, 3, 0, 1, 0, 0], 'T': [1, 5, 0, 0, 0, 1, 1, 6]}
    assert expected == profile_matrix


class TestUnite:
    def test_unite_folder(self, base_data_path):
        unite_fasta_in_folder(base_data_path + "unite", file_name_part_to_id=".")
        actual = read_fasta(base_data_path + "unite.fasta")
        assert actual == [('1.id1', 'aaa'), ('1.id2', 'bbb'), ('2.id2__1', 'ccc'), ('3.id1__1', 'abc'), ('3.id5', 'ggg')]
