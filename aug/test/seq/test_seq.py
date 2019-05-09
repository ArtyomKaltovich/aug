import textwrap

import pytest

from aug.comb.comb import gen_substrings
from aug.heredity.Phenotype import *
from aug.heredity.heredity import n_expected_dominant_phenotype
from aug.seq.seq import *
from aug.test import FLOAT_EQUALITY_ACCURACY


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


def test_all_possible_gene_transcription():
    string = "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
    expected = {"MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"}
    actual = all_possible_gene_transcription(string)
    assert expected == actual


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


def test_failure_array():
    dna = "ACACAC"
    expected = [0, 0, 1, 2, 3, 4]
    assert expected == failure_array(dna)
    dna = "ACGTACGACGTATT"
    expected = [0, 0, 0, 0, 1, 2, 3, 1, 2, 3, 4, 5, 0, 0]
    assert expected == failure_array(dna)


def test_dominant_phenotype_prob():
    p = PhenotypeHeredityTable()
    assert 0.75 == p.dominant_phenotype_prob(("Aa", "Aa"))
    assert 0.0 == p.dominant_phenotype_prob(("aa", "aa"))
    assert 0.5 == p.dominant_phenotype_prob(("Aa", "aa"))
    assert 1.0 == p.dominant_phenotype_prob(("AA", "AA"))


def test_n_expected_dominant_phenotype():
    assert 0 == n_expected_dominant_phenotype([0, 0, 0, 0, 0, 1])
    assert 5 == n_expected_dominant_phenotype([1, 0, 0, 0, 0, 0], n_children=5)
    assert 7.5 == n_expected_dominant_phenotype([0, 0, 0, 10, 0, 0])


def test_gen_substrings():
    actual = list(gen_substrings("DNA", 3))
    expected = [ "D", "DD", "DDD", "DDN", "DDA", "DN", "DND", "DNN", "DNA", "DA", "DAD", "DAN", "DAA",
                 "N", "ND", "NDD", "NDN", "NDA", "NN", "NND", "NNN", "NNA", "NA", "NAD", "NAN", "NAA",
                 "A", "AD", "ADD", "ADN", "ADA", "AN", "AND", "ANN", "ANA", "AA", "AAD", "AAN", "AAA"]
    assert actual == expected


def test_adjacency_list():
    expected = [('Rosalind_0498', 'Rosalind_2391'), ('Rosalind_0498', 'Rosalind_0442'), ('Rosalind_2391', 'Rosalind_2323')]
    actual = adjacency_list("test_data_files/test_adjacency_list.txt", k=3)
    assert len(expected) == len(actual)
    for val in expected:
        assert val in actual


def test_find_protein_motif_by_shorthand():
    protein = """
        MLGVLVLGALALAGLGFPAPAEPQPGGSQCVEHDCFALYPGPATFLNASQICDGLRGHLM
        TVRSSVAADVISLLLNGDGGVGRRRLWIGLQLPPGCGDPKRLGPLRGFQWVTGDNNTSYS
        RWARLDLNGAPLCGPLCVAVSAAEATVPSEPIWEEQQCEVKADGFLCEFHFPATCRPLAV
        EPGAAAAAVSITYGTPFAARGADFQALPVGSSAAVAPLGLQLMCTAPPGAVQGHWAREAP
        GAWDCSVENGGCEHACNAIPGAPRCQCPAGAALQADGRSCTASATQSCNDLCEHFCVPNP
        DQPGSYSCMCETGYRLAADQHRCEDVDDCILEPSPCPQRCVNTQGGFECHCYPNYDLVDG
        ECVEPVDPCFRANCEYQCQPLNQTSYLCVCAEGFAPIPHEPHRCQMFCNQTACPADCDPN
        TQASCECPEGYILDDGFICTDIDECENGGFCSGVCHNLPGTFECICGPDSALARHIGTDC
        DSGKVDGGDSGSGEPPPSPTPGSTLTPPAVGLVHSGLLIGISIASLCLVVALLALLCHLR
        KKQGAARAKMEYKCAAPSKEVVLQHVRTERTPQRL"""
    protein = textwrap.dedent(protein).replace("\n", "")
    assert [46, 114, 115, 381, 408] == find_protein_motif_by_shorthand(protein, "N{P}[ST]{P}")


def test_transition_transversion_ratio():
    dna1 = "GCAACGCACAACGAAAACCCTTAGGGACTGGATTATTTCGTGATCGTTGTAGTTATTGGAAGTACGGGCATCAACCCAGTT"
    dna2 = "TTATCTGACAAAGAAAGCCGTCAACGGCTGGATAATTTCGCGATCGTGCTGGTTACTGGCGGTACGAGTGTTCCTTTGGGT"
    assert 1.214 == pytest.approx(transition_transversion_ratio(dna1, dna2), FLOAT_EQUALITY_ACCURACY)
