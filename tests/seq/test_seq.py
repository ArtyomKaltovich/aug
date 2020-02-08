import random
import textwrap

import pytest

from aug.comb.comb import gen_substrings
from aug.heredity.Phenotype import *
from aug.heredity.heredity import n_expected_dominant_phenotype
from aug.seq.seq import *
from tests import FLOAT_EQUALITY_ACCURACY
from tests.utils import random_string


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
    assert expected == pytest.approx(gc_rate(param, percent=True))


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
    assert 821.392 == pytest.approx(calculate_protein_mass("SKADYEK"))


def test_dominant_probability():
    assert 0.78333 == pytest.approx(dominant_probability(2, 2, 2))


def test_consensus():
    data = ["ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC", "TTGGAACT", "ATGCCATT", "ATGGCACT", ]
    assert "ATGCAACT" == consensus(data)


def test_consensus_precalculated():
    data = {'A': [5, 1, 0, 0, 5, 5, 0, 0], 'C': [0, 0, 1, 4, 2, 0, 6, 1], 'G': [1, 1, 6, 3, 0, 1, 0, 0],
            'T': [1, 5, 0, 0, 0, 1, 1, 6]}
    assert "ATGCAACT" == consensus(precalculated_profile=data)


def test_n_reverse_translation():
    assert 12 == n_reverse_translation("MA")


def test_find_reverse_palindromes():
    expected = [(4, 6), (5, 4), (6, 6), (7, 4), (17, 4), (18, 4), (20, 6), (21, 4)]
    assert expected == find_reverse_palindromes("TCAATGCATGCGGGTCTATATGCAT", zero_based=False)


def test_dna_probability():
    assert 1.831e-06 == pytest.approx(dna_probability("ACGATACAA", 0.129))
    assert -5.737 == pytest.approx(dna_probability("ACGATACAA", 0.129, return_log=True))


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


@pytest.mark.parametrize("method, score_coef, score_delta",
                         [[alignments.Levinshtein(), 1, 0],
                          [alignments.NeedlemanWunsch(), -1, 0],
                          [alignments.NeedlemanWunsch(gap_start=10), -1, 10],
                          [alignments.NeedlemanWunsch(gap_start=-10), -1, -10]])
def test_align_empty(random_seed, method, score_coef, score_delta):
    string = random_string(min_len=1)
    score = score_coef * len(string) + score_delta
    assert align(string, "", reconstruct_answer=True, method=method, swap_case_on_mismatch=False) == \
                ((string, "-" * len(string)), score)
    assert align("", string, reconstruct_answer=True, method=method, swap_case_on_mismatch=False) == \
                (("-" * len(string), string), score)


@pytest.mark.parametrize("method",
                         [alignments.Levinshtein(),
                          alignments.NeedlemanWunsch(),
                          alignments.NeedlemanWunsch(gap_start=10),
                          alignments.NeedlemanWunsch(gap_start=-10)])
def test_align_random_strings(random_seed, method):
    # seq part between two gaps should be presented in the seq.
    string1 = random_string(min_len=5)
    string2 = random_string(min_len=5)
    (alignment1, alignment2), score = align(string1, string2, reconstruct_answer=True, method=method,
                                            swap_case_on_mismatch=False)
    substring1 = random.choice([a for a in alignment1.split("-") if a])
    assert substring1 in string1
    substring2 = random.choice([a for a in alignment2.split("-") if a])
    assert substring2 in string2


def test_edit_distance():
    assert 0 == edit_distance("ab", "ab")
    assert 2 == edit_distance("", "ab")
    assert 3 == edit_distance("short", "ports")
    assert 5 == edit_distance("editing", "distance")


@pytest.mark.parametrize("seq1, seq2, alignment1, alignment2, score",
                         [["ab", "ab", "ab", "ab", 0],
                          ["editing", "distance", "edi-tin-g", "-distance", 5],
                          ["banana", "ananas", "banana-", "-ananas", 2],
                          ["ananas", "banana", "-ananas", "banana-", 2],
                          ["ananas", "anna", "ananas", "an-na-", 2],
                          ["bamboo", "baobab", "bamboo", "baobab", 3],
                          ["zzzbamboo", "baobab", "zzzbamboo", "---baobab", 6],
                          ])
def test_edit_distance_restore_answer(seq1, seq2, alignment1, alignment2, score):
    assert edit_distance(seq1, seq2, reconstruct_answer=True, swap_case_on_mismatch=False) \
           == ((alignment1, alignment2), score)


def test_align_needleman_wunsch_as_edit_distance(random_seed):
    # tests that Needleman Wunsch works the same as edit distance with the same scores.
    string1 = random_string()
    string2 = random_string()
    (expected1, expected2), _ = edit_distance(string1, string2, reconstruct_answer=True,
                                              swap_case_on_mismatch=False)
    method = alignments.NeedlemanWunsch(match_score=0, mismatch_score=-1, gap_score=-1)
    (actual1, actual2), _ = align(string1, string2, reconstruct_answer=True, method=method, swap_case_on_mismatch=False)
    assert (expected1, expected2) == (actual1, actual2)


@pytest.mark.parametrize("seq1, seq2",
                         [["axc", "aabcc"],
                          ["MGTSADNALAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA",
                           "MLKREVLRDRKVFGNPIACRQEVFRWCMRYNTHRRHSWCNLVAPDVFETETSATLTKAT"],
                          ["MLQTPCPQSGVRVPLAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA",
                           "MKREVLQDAACWPDEATCRRQVFRWAVRYNTRRRHSWCGYLSPSTYEARWAATLPTAA"],
                          ["MKREVLQDAACWPDEATCRRQVFRWAVRYNTRRRHSWCGYLSPSTYEARWAATLPTAA",
                           "MLKREVLRDRKVFGNPIACRQEVFRWCMRYNTHRRHSWCNLVAPDVFETETSATLTKAT"]])
def test_align_affine_gap(seq1, seq2):
    #  with big penalty for opening gap les gap should be presented
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=1)
    (more_gaps1, more_gaps2), _ = align(seq1, seq2, reconstruct_answer=True, method=method)
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=-10)
    (less_gaps1, less_gaps2), _ = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert len(more_gaps1.split("-")) >= len(less_gaps1.split("-"))
    assert len(more_gaps2.split("-")) >= len(less_gaps2.split("-"))


def test_needleman_wunsh_by_hamming_distance_random_strings(random_seed):
    seq1 = random_string()
    seq2 = random_string()
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method, swap_case_on_mismatch=False)
    assert score == len(line1) - 2 * hamming_distance(line1, line2), (seq1, seq2)


def def_test_align_with_distance_matrix():
    seq1 = "XXXYZ"
    seq2 = "ЙYXXX"

    dist_matrix = {'X': {'Й': -5, 'Z': -5, 'Y': -5, 'X': 5},
                   'Y': {'Й': -5, 'Z': -5, 'Y': 5, 'X': -5},
                   'Z': {'Й': -15, 'Z': 8, 'Y': -5, 'X': -5},
                   'Й': {'Й': 8, 'Z': -15, 'Y': -5, 'X': -5}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix, gap_score=-10)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert ((line1, line2), score) == (("XXXYZ", "йyXxx"), -15)

    dist_matrix = {'X': {'Й': -5, 'Z': -5, 'Y': -5, 'X': 15},
                   'Y': {'Й': -5, 'Z': -5, 'Y': -5, 'X': -5},
                   'Z': {'Й': -15, 'Z': 8, 'Y': -5, 'X': -5},
                   'Й': {'Й': 8, 'Z': -15, 'Y': -5, 'X': -5}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix, gap_score=-10)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert ((line1, line2), score) == (("--XXXYZ", "йyXXX--"), 5)


def test_align_with_gap_in_distance_matrix():
    seq1 = "CBBBC"
    seq2 = "ABBBA"

    dist_matrix = {'A': {'A': 5, 'B': -5, 'C': -5, '-': -5},
                   'B': {'A': -5, 'B': 5, 'C': -15, '-': -5},
                   'C': {'A': -5, 'B': -15, 'C': 5, '-': -5},
                   '-': {'A': -5, 'B': -5, 'C': 3, '-': 0}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert ((line1, line2), score) == (("-CBBB-C", "a-BBBa-"), 11)

    dist_matrix = {'A': {'A': 5, 'B': -5, 'C': -5, '-': -5},
                   'B': {'A': -5, 'B': -15, 'C': -15, '-': -5},
                   'C': {'A': -5, 'B': -15, 'C': 5, '-': -5},
                   '-': {'A': -5, 'B': -5, 'C': 3, '-': 0}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert ((line1, line2), score) == (("CB---BBC", "-abbb-a-"), -24)


def test_random_align_with_affine_gap_in_distance_matrix(random_seed):
    # if we specify positive gap start score there should be more gaps in alignment
    seq1 = random_string(min_len=3, alphabet=["A", "B", "C"])
    seq2 = random_string(min_len=3, alphabet=["A", "B", "C"])

    dist_matrix = {'A': {'A': 5, 'B': -5, 'C': -5, '-': -5},
                   'B': {'A': -5, 'B': 5, 'C': -15, '-': -5},
                   'C': {'A': -5, 'B': -15, 'C': 5, '-': -5},
                   '-': {'A': -5, 'B': -5, 'C': 3, '-': 0}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix)
    (alignment1, alignment2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix, gap_start=20)
    (gap_start_alignment1, gap_start_alignment2), gap_start_score \
        = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert len(alignment1.split("-")) <= len(gap_start_alignment1.split("-"))
    assert len(alignment2.split("-")) <= len(gap_start_alignment2.split("-"))


@pytest.mark.parametrize("name, value", [
    ["match_score", 10],
    ["mismatch_score", 10]
])
def test_score_matrix_and_mis_match_score_passed_as_param(name, value):
    kwargs = {name: value}
    with pytest.raises(ValueError):
        alignments.NeedlemanWunsch(score_matrix=blosum62, gap_score=-10, **kwargs)


def test_needleman_wunsh_affine_gap():
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=-10)
    (line1, line2), score = align("axc", "aabcc", reconstruct_answer=True, method=method)
    assert (('a--xc', 'aABCc'), -11) == ((line1, line2), score)
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=1)
    (line1, line2), score = align("axc", "aabcc", reconstruct_answer=True, method=method)
    assert (('-a-x-c', 'AaB-Cc'), 2) == ((line1, line2), score)
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=2)
    (line1, line2), score = align("axc", "aabcc", reconstruct_answer=True, method=method)
    assert (('-a-x-c-', 'A-A-BcC'), 7) == ((line1, line2), score)


@pytest.mark.parametrize("seq1, seq2, alignment1, alignment2, score",
                         [["ACC", "AACCC", "  ACC ", "A ACC C", 3],
                          ["TGTTACGG", "GGTTGACTA", "  GTT-AC GG", "G GTTgAC TA", 4]])
def test_local_alignment(seq1, seq2, alignment1, alignment2, score):
    method = alignments.SmithWaterman(match_score=1, mismatch_score=-1, gap_score=-1)
    (line1, line2), actual_score = align(seq1, seq2, reconstruct_answer=True, method=method)
    assert ((line1, line2), actual_score) == ((alignment1, alignment2), score)


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
    expected = ["D", "DD", "DDD", "DDN", "DDA", "DN", "DND", "DNN", "DNA", "DA", "DAD", "DAN", "DAA",
                "N", "ND", "NDD", "NDN", "NDA", "NN", "NND", "NNN", "NNA", "NA", "NAD", "NAN", "NAA",
                "A", "AD", "ADD", "ADN", "ADA", "AN", "AND", "ANN", "ANA", "AA", "AAD", "AAN", "AAA"]
    assert actual == expected


def test_adjacency_list(base_data_path):
    expected = [('Rosalind_0498', 'Rosalind_2391'), ('Rosalind_0498', 'Rosalind_0442'),
                ('Rosalind_2391', 'Rosalind_2323')]
    actual = adjacency_list(base_data_path + r"test_adjacency_list.txt", k=3)
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
    assert 1.214 == pytest.approx(transition_transversion_ratio(dna1, dna2))


#def test_longest_common_substring():
#    assert ['ba', 'ab'] == longest_common_substring(["abba", "baba"])
#    assert ['anana'] == longest_common_substring(["banana", "ananas"])
#    assert ['ban'] == longest_common_substring(["banana", "ban"])
#
#
#def test_longest_common_substring3():
#    assert ['ba'] == longest_common_substring(["baobab", "bamboo", "banana"])
#    assert ['an'] == longest_common_substring(["ananas", "ban", "banana"])\
#           == longest_common_substring(["banana", "ananas", "ban"])
#    assert 'AC' in longest_common_substring(["GATTACA", "TAGACCA", "ATACA"])
