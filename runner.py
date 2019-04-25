from utils import *
from Sets import Sets
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
    reader = read_file()
    n = int(next(reader))
    sets = Sets(n)
    with open("data/answer.txt", "w") as output:
        for line in reader:
            a, b = map(int, line.split())
            sets.unite(a - 1, b - 1)
        print(sets.n_disjoint - 1, file=output)
