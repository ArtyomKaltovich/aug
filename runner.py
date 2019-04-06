from utils import *

def read_digits(line):
    i, j = map(int, line.split())
    return i, j


def for_fasta():
    for id, string in fasta_file_iter("data/rosalind_dna.txt"):
        result = find_reverse_palindromes(string, zero_based=False)
    for r in result:
        print(*r)


def read_file():
    with open("data/rosalind_dna.txt", "r") as f:
        line = f.readline().strip()
        #i, j = read_digits(dna)

        print(find_reverse_palindromes(line))


if __name__ == '__main__':
    gene, *intrones = fasta_file_iter("data/rosalind_dna.txt")
    intrones = [string for id, string in intrones]
    print(gene_to_protein(gene[1], intrones))