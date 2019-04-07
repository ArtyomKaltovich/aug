from utils import *

def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        i = list(map(int, line.split()))
    return i


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
    n, *_ = read_digits()
    permutations = list(signed_permutation(n))
    print(len(permutations))
    for p in permutations:
        print(*p)
