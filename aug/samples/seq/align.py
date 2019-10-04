import random
from collections import defaultdict

from aug.seq import alignments
from aug.seq.seq import edit_distance, align, hamming_distance


def simple_needleman_wunsch_sample():
    seq1 = "MGTSADNALAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA"
    seq2 = "MLQTPCPQSGVRVPLAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPD"
    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-0.499)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)


def affine_gap_needleman_wunsch_sample():
    seq1 = "TCCCAGTTATGTCAGGGGACACGAGCATGCAGAGAC"
    seq2 = "AATTGCCGCCGTCGTTTTCAGCAGTTATGTCAGATC"

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_start=0, gap_score=-1)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_start=-100, gap_score=-0.01)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_start=0.5, gap_score=-0.3)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)


def distance_matrix_sample():
    seq1 = "XXXYZ"
    seq2 = "ЙYXXX"

    dist_matrix = {'X': {'-': -3, 'Й': -5, 'Z': -5, 'Y': -5, 'X': 5},
                   'Y': {'-': -5, 'Й': -5, 'Z': -5, 'Y': 5, 'X': -5},
                   'Z': {'-': 3, 'Й': -15, 'Z': 8, 'Y': -5, 'X': -5},
                   'Й': {'-': -3, 'Й': 8, 'Z': -15, 'Y': -5, 'X': -5},
                   '-': {'-': 0, 'Й': -3, 'Z': 3, 'Y': -5, 'X': -3}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    dist_matrix = {'X': {'-': 33, 'Й': -5, 'Z': -5, 'Y': -5, 'X': 5},
                   'Y': {'-': -5, 'Й': -5, 'Z': -5, 'Y': 5, 'X': -5},
                   'Z': {'-': 3, 'Й': -15, 'Z': 8, 'Y': -5, 'X': -5},
                   'Й': {'-': -3, 'Й': 8, 'Z': -15, 'Y': -5, 'X': -5},
                   '-': {'-': 0, 'Й': -3, 'Z': 3, 'Y': -5, 'X': 33}}

    method = alignments.NeedlemanWunsch(score_matrix=dist_matrix)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)


def local_alignment_sample():
    seq1 = "ACC"
    seq2 = "AACCC"
    #seq1 = "ALAESFNSALKQTPCPQSGVRVPLREVLQDRKVFDNHLVCYEAAA"
    #seq2 = "MLQTPCPQSGVRVPLAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA"  # https://www.uniprot.org/uniprot/W5Y176.fasta

    method = alignments.NeedlemanWunsch(match_score=3, mismatch_score=-3, gap_score=-2)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    method = alignments.SmithWaterman(match_score=1, mismatch_score=-1, gap_score=-1)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    seq1 = "TGTTACGG"
    seq2 = "GGTTGACTA"
    method = alignments.SmithWaterman(match_score=1, mismatch_score=-1, gap_score=-1)
    (line1, line2), score = align(seq1, seq2, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)




def test():
    seq1 = "MGTSADNALAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA"  # https://www.uniprot.org/uniprot/W5Y845.fasta
    seq2 = "MLQTPCPQSGVRVPLAESFNSALKREVLQDRKVFDNHLVCRREVFYWCTRYNTHRLHTWCGYLSPDDYEAAA"  # https://www.uniprot.org/uniprot/W5Y176.fasta
    seq3 = "MKREVLQDAACWPDEATCRRQVFRWAVRYNTRRRHSWCGYLSPSTYEARWAATLPTAA"  # https://www.uniprot.org/uniprot/A0A345NSW8.fasta
    seq4 = "MLKREVLRDRKVFGNPIACRQEVFRWCMRYNTHRRHSWCNLVAPDVFETETSATLTKAT" # https://www.uniprot.org/uniprot/A0A127NSH0.fasta

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=-1)
    (line1, line2), score = align(seq4, seq3, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    method = alignments.NeedlemanWunsch(match_score=1, mismatch_score=-1, gap_score=-1, gap_start=-10)
    (line1, line2), score = align(seq4, seq3, reconstruct_answer=True, method=method)
    print(line1)
    print(line2)
    print(score)

    (line1, line2), score = edit_distance(seq3, seq4, reconstruct_answer=True)
    print(line1)
    print(line2)
    print(score)


if __name__ == '__main__':
    #simple_needleman_wunsch_sample()
    #affine_gap_needleman_wunsch_sample()
    distance_matrix_sample()
    #local_alignment_sample()
