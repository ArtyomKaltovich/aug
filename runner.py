from utils import *
from Sets import Sets
from time import perf_counter
import math
from scipy.special import comb

def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        n, k = list(map(int, line.split()))
    return n, k


def for_fasta():
    for id, string in fasta_file_iter("data/rosalind_data.txt"):
        result = failure_array(string)
    #for r in result:
    print(*result)


def read_file():
    with open("data/rosalind_data.txt", "r") as f:
        for line in f:
            yield line.strip()
        #i, j = read_digits(dna)


if __name__ == '__main__':
    #file=open("data/answer.txt", "w")
    #data = read_fasta("data/rosalind_data.txt", without_id=True)

