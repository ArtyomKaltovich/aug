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


@pytest.mark.skip("file operations")
def test_fasta_file_do():
    def helper(id, string):
        return id, string

    actual = fasta_data_do(r"data/test_fasta.txt", helper, [])
    expected = [("Rosalind_6404", "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("Rosalind_5959", "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("Rosalind_0808", "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == actual


@pytest.mark.skip("file operations")
def test_read_fasta():
    expected = [("Rosalind_6404", "CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCCTCCCACTAATAATTCTGAGG"),
                ("Rosalind_5959", "CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCTATATCCATTTGTCAGCAGACACGC"),
                ("Rosalind_0808", "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT")]
    assert expected == read_fasta(r"data/test_fasta.txt")


@pytest.mark.skip("file operations")
def test_fasta_file_do_max():
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
    assert 60.919540 == pytest.approx(gc_max, FLOAT_EQUALITY_ACCURACY)


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
    assert [6] == find_motif("AGGCTGAATTCAGATCAC", "GAATTC", zero_based=False)


def test_find_motif_long():
    assert [237] == find_motif(
        """AAGGCAACATTTCTTAGTTATATATGCTTGTAGTGAAGAAAGATGTGAAAGTCTGACAAGAGAACAAGAC
            GAAGGAGGAGTCTTTCTCCAAGTCTTCAACATTGCAGAATCTGATGCATATGAACCCATTTTCTCTACAA
            AATGTTGCAACCCTAGAGAGCAAAACAAAACATACCCATAATCAGAAATGATCTGACGAAAATCGAGTTA
            CAATACACAAGAGAACATTTTTTTTAGAATTCTCAGATATTAAAAATGACACAGAAAGCTTTATGCTTTT
            TCCTCTTAAAAGACTAAACAAGTTGAAATCTAGAGAAAGAACTGACCAACCTGAGACAACGAGAGAGACT
            TGAGAGATTTCTTCGGCACTTACTATTAGATCTAGGGTTTAGATACCATTTATATAGAGAAAGTTTTAGA
            GTTGCACAAAACATAAATTAATGTGTTAGAATGGGCCTAAAGCTACAAAGCTGGCCTGGTTTTGTTTTAA
            ATTGTTGGTTTCATGGACATTTTCGACATCTTCGAACATGTTATTTTTTGAGACTATGCAAACTTGAAGC
            TCTTTACTCGAGTTGAAATCGTATGACTTATAGTGAAATTGTACATTTGGTTTCGATTTTTCTTTTACAC
            TCTTTCTTCTTTGAGCCGGTAAATTTGGAATTTTTCTTCATAGTGGAATCATATGCTGTTTTTTTTTTTT
            ATAGTAAACGTTACAAGAATGAATGGTAACTTTATCCAAAAAAAAAGAATCATATTATTTTGAAATGATT
            TTAAGTAAATTCTAGGTTCAATAACATAAGATTTGAGACTAAATTTAAAATTTCTTAGTAAAATATATGA
            TTTTTTTATAAATACCTATAAAATTAGTAATTAACAATACGGATTACGTACTGAATCAAACCCTTTGTAT
            TTTGTTTTTCCTAGAAATAAGTGTAGATTTTTGGAATTTTGCATTAATTAATCACTTCTTGGGTCTGAAA
            GGCTAAAACAAAAGGAACCGAAAGAGAATGTTCTCTCTGTCTTTATCTTCCACTTCCACTTCCAGGTCGC
            GTTGCTTCACTCTCCATTGCAAAGAGAGGTCTCTGCGATTTCTGCAACTCACCCCTGAAACCTTCTTAAT
            TTACTTCAACTGCCGCTATACCTAAAAACTTCATCTTTCTCCTCTGAGCTATG""".replace("\n", "").replace(" ", ""),
        motif="GAATTC", zero_based=False)


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


@pytest.mark.skip("file operations")
def test_adjacency_list():
    expected = [('Rosalind_0498', 'Rosalind_2391'), ('Rosalind_0498', 'Rosalind_0442'), ('Rosalind_2391', 'Rosalind_2323')]
    actual = adjacency_list("data/test_adjacency_list.txt", k=3)
    assert len(expected) == len(actual)
    for val in expected:
        assert val in actual


def test_dna_probability():
    assert 1.831e-06 == pytest.approx(dna_probability("ACGATACAA", 0.129), FLOAT_EQUALITY_ACCURACY)
    assert -5.737 == pytest.approx(dna_probability("ACGATACAA", 0.129, return_log=True), FLOAT_EQUALITY_ACCURACY)


def test_find_spliced_motif():
    dna = "ACGTACGTGACG"
    motif = "GTA"
    actual = find_spliced_motif(dna, motif, zero_based=False)
    for a, i in zip(actual, range(len(motif))):
        assert dna[a - 1] == motif[i]


def test_find_spliced_motif_non_found():
    dna = "ACGTACGTGACG"
    motif = "123"
    assert -1 == find_spliced_motif(dna, motif, zero_based=False)


def test_edit_distance():
    assert 0 == edit_distance("ab", "ab")
    assert 2 == edit_distance("", "ab")
    assert 3 == edit_distance("short", "ports")
    assert 5 == edit_distance("editing", "distance")


def test_enumerate_kmers():
    assert "ACGT" == "".join(enumerate_kmers("ACGT", 1))
    assert "AAACCACC" == "".join(enumerate_kmers("AC", 2))


def test_string_to_kmers():
    assert [] == list(string_to_kmers("", 2))
    assert ["aaaa", "bbbb", "cccc"] == list(string_to_kmers("aaaabbbbcccc", 4))


def test_kmers_composition():
    assert [0] * 16 == list(kmers_composition("", 2))
    assert [2, 1, 1, 2] == list(kmers_composition("ACGATT", 1))


def test_distance_matrix():
    dnas = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
    assert [[0, 0.4, 0.1, 0.1], [0.4, 0, 0.4, 0.3], [0.1, 0.4, 0, 0.2], [0.1, 0.3, 0.2, 0]] == distance_matrix(dnas)


def test_distance_matrix_non_relative():
    dnas = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
    assert [[0, 4, 1, 1], [4, 0, 4, 3], [1, 4, 0, 2], [1, 3, 2, 0]] == distance_matrix(dnas, relative=False)


def test_distance_matrix_ndarray():
    dnas = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
    expected = np.asarray([[0, 0.4, 0.1, 0.1], [0.4, 0, 0.4, 0.3], [0.1, 0.4, 0, 0.2], [0.1, 0.3, 0.2, 0]])
    np.testing.assert_array_equal(expected, distance_matrix(dnas, as_ndarray=True))
