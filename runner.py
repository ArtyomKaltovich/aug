from aug.comb.comb import gen_substrings
from aug.data.fasta import *
from aug.data.net import request_protein_sequence
from aug.seq.seq import *


def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        n = list(map(int, line.split()))
    return n


def for_fasta():
    for id, string in fasta_file_iter("data/rosalind_data.txt"):
        result = all_possible_gene_transcription(string)
    for r in result:
        print(r, file=open("data/answer.txt", "a"))


def read_file():
    with open("data/rosalind_data.txt", "r") as f:
        for line in f:
            yield line.strip()
        #i, j = read_digits(dna)


if __name__ == '__main__':
    #file=open("data/answer.txt", "w")
    data = read_fasta("data/rosalind_data.txt", without_id=True)
    print(transition_transversion_ratio(*data))

