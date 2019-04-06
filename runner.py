from utils import *

def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        i, j = map(int, line.split())
    return i, j


def for_fasta():
    for id, string in fasta_file_iter("data/rosalind_data.txt"):
        result = find_reverse_palindromes(string, zero_based=False)
    for r in result:
        print(*r)


def read_file():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        #i, j = read_digits(dna)

        print(find_reverse_palindromes(line))


if __name__ == '__main__':
    j, i = read_digits()
    print(j, i)
    result = independent_alleles(i, j)
    print(round(result, 3))