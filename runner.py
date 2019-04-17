from utils import *
from time import perf_counter
import math

def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        n, k = list(map(int, line.split()))
    return n, k


def for_fasta():
    for id, string in fasta_file_iter("data/rosalind_data.txt"):
        result = find_reverse_palindromes(string, zero_based=False)
    for r in result:
        print(*r)


def read_file():
    with open("data/rosalind_data.txt", "r") as f:
        for line in f:
            yield line.strip()
        #i, j = read_digits(dna)


if __name__ == '__main__':
    for id, string in fasta_file_iter("data/rosalind_data.txt"):
        print(*count_kmers(string, 4), file=open("data/answer.txt", "w"))
