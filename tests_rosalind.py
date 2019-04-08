import pytest
from utils import *

FLOAT_EQUALITY_ACCURACY = 0.001


def test_dna_to_rna():
    param = "GATGGAACTTGACTACGTAAATT"
    expected = "GAUGGAACUUGACUACGUAAAUU"
    assert expected == dna_to_rna(param)


def test_rna_to_protein():
    assert "MAMAPRTEINSTRING" == rna_to_protein("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")
    assert "MATE" == rna_to_protein("AGAUGGCCACCGAGUAA", start=start_with_start_codon)


def test_rna_to_protein_only_stop():
    assert "" == rna_to_protein("UAA")


def test_gene_to_protein():
    gene = "ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"
    intrones = "ATCGGTCGAA", "ATCGGTCGAGCGTGT"
    assert "MVYIADKQHVASREAYGHMFKVCA" == gene_to_protein(gene, intrones)


def test_reverse_complement():
    param = "AAAACCCGGT"
    expected = "ACCGGGTTTT"
    assert expected == reverse_complement(param)


@pytest.mark.skip
def test_fasta_file_do():
    def helper(id, string):
        return id, string

    actual = fasta_data_do(r"data/test_fasta.txt", helper, [])
    expected = [("Rosalind_6404", "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("Rosalind_5959", "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("Rosalind_0808", "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == actual


@pytest.mark.skip
def test_fasta_file_do_sort():
    id_max = None
    gc_max = None

    def helper(id, string):
        nonlocal id_max, gc_max
        gc = gc_rate(string, procent=True)
        if not gc_max or gc > gc_max:
            id_max = id
            gc_max = gc

    fasta_data_do(r"data/test_fasta.txt", helper)
    assert "Rosalind_0808" == id_max
    assert 60.919540 == pytest.approx(gc_max)


def test_gc_rate():
    param = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"
    expected = 60.9195
    assert expected == pytest.approx(gc_rate(param, procent=True), FLOAT_EQUALITY_ACCURACY)


def test_hamming_distance():
    assert 7 == hamming_distance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")


def test_find_motif():
    assert [1, 3, 9] == find_motif("GATATATGCATATACTT", "ATAT", zero_based=True)
    assert [2, 4, 10] == find_motif("GATATATGCATATACTT", "ATAT", zero_based=False)
    assert [0, 10, 17] == find_motif("find_me789find_mefind_mefind_m", "find_me")
    assert [3] == find_motif("012find_me", "find_me")


def test_rabbits_recurrence():
    assert 12586269025 == rabbits_recurrence(50)
    assert 19 == rabbits_recurrence(5, 3)


def test_dying_rabbits():
    assert 12586269025 == dying_rabbits(50, 51)
    assert 1 == dying_rabbits(1, 3)
    assert 1 == dying_rabbits(2, 3)
    assert 2 == dying_rabbits(3, 3)
    assert 2 == dying_rabbits(4, 3)
    assert 3 == dying_rabbits(5, 3)
    assert 4 == dying_rabbits(6, 3)
    assert 5 == dying_rabbits(7, 3)


def test_calculate_protein_mass():
    assert 821.392 == pytest.approx(calculate_protein_mass("SKADYEK"), FLOAT_EQUALITY_ACCURACY)


def test_dominant_probability():
    assert 0.78333 == pytest.approx(dominant_probability(2, 2, 2), FLOAT_EQUALITY_ACCURACY)


def test_fasta_profile():
    profile_matrix = None
    for id, dna in fasta_file_iter("data/test_profile_fasta.txt"):
        profile_matrix = profile(dna, profile_matrix)
    expected = {'A': [5, 1, 0, 0, 5, 5, 0, 0], 'C': [0, 0, 1, 4, 2, 0, 6, 1], 'G': [1, 1, 6, 3, 0, 1, 0, 0], 'T': [1, 5, 0, 0, 0, 1, 1, 6]}
    assert expected == profile_matrix


def test_consensus():
    data = ["ATCCAGCT","GGGCAACT","ATGGATCT","AAGCAACC","TTGGAACT","ATGCCATT","ATGGCACT",]
    assert "ATGCAACT" == consensus(data)


def test_consensus_precalculated():
    data = {'A': [5, 1, 0, 0, 5, 5, 0, 0], 'C': [0, 0, 1, 4, 2, 0, 6, 1], 'G': [1, 1, 6, 3, 0, 1, 0, 0], 'T': [1, 5, 0, 0, 0, 1, 1, 6]}
    assert "ATGCAACT" == consensus(precalculated_profile=data)


def test_n_reverse_translation():
    assert 12 == n_reverse_translation("MA")


def test_find_reverse_palindromes():
    expected = [(4, 6), (5, 4), (6, 6), (7, 4), (17, 4), (18, 4), (20, 6), (21, 4)]
    assert expected == find_reverse_palindromes("TCAATGCATGCGGGTCTATATGCAT", zero_based=False)
