from collections import Counter
from functools import lru_cache
from typing import List, Union, Dict, Tuple, Collection
from scipy.special import comb
from itertools import permutations, combinations_with_replacement, product

complement_map = {"A": "T", "C": "G", "G": "C", "T": "A"}

STOP_CODON = object()  # dummy object for stop codon
rna_codon_table = { "UUU": "F",         "CUU": "L",     "AUU": "I",      "GUU": "V",
                    "UUC": "F",         "CUC": "L",     "AUC": "I",      "GUC": "V",
                    "UUA": "L",         "CUA": "L",     "AUA": "I",      "GUA": "V",
                    "UUG": "L",         "CUG": "L",     "AUG": "M",      "GUG": "V",
                    "UCU": "S",         "CCU": "P",     "ACU": "T",      "GCU": "A",
                    "UCC": "S",         "CCC": "P",     "ACC": "T",      "GCC": "A",
                    "UCA": "S",         "CCA": "P",     "ACA": "T",      "GCA": "A",
                    "UCG": "S",         "CCG": "P",     "ACG": "T",      "GCG": "A",
                    "UAU": "Y",         "CAU": "H",     "AAU": "N",      "GAU": "D",
                    "UAC": "Y",         "CAC": "H",     "AAC": "N",      "GAC": "D",
                    "UAA": STOP_CODON,  "CAA": "Q",     "AAA": "K",      "GAA": "E",
                    "UAG": STOP_CODON,  "CAG": "Q",     "AAG": "K",      "GAG": "E",
                    "UGU": "C",         "CGU": "R",     "AGU": "S",      "GGU": "G",
                    "UGC": "C",         "CGC": "R",     "AGC": "S",      "GGC": "G",
                    "UGA": STOP_CODON,  "CGA": "R",     "AGA": "R",      "GGA": "G",
                    "UGG": "W",         "CGG": "R",     "AGG": "R",      "GGG": "G",
}
protein_n_codons_table = Counter(rna_codon_table.values())

monoisotopic_mass_table = {
    "A":  71.03711,
    "C":  103.00919,
    "D":  115.02694,
    "E":  129.04259,
    "F":  147.06841,
    "G":  57.02146,
    "H":  137.05891,
    "I":  113.08406,
    "K":  128.09496,
    "L":  113.08406,
    "M":  131.04049,
    "N":  114.04293,
    "P":  97.05276,
    "Q":  128.05858,
    "R":  156.10111,
    "S":  87.03203,
    "T":  101.04768,
    "V":  99.06841,
    "W":  186.07931,
    "Y":  163.06333,
}


def dna_to_rna(dna: str):
    """An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

    Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all
        occurrences of 'T' in t with 'U' in u.
    >>> dna_to_rna("GCAT")
    'GCAU'
    """
    return dna.replace("T", "U")


def rna_to_protein(rna: str, to_str=True):
    """The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
        Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
        Return: The protein string encoded by s.
    """
    result = [""] * (len(rna) // 3)
    for i in range(len(result)):
        elem = rna_codon_table[rna[3 * i: 3 * i + 3]]
        if elem != STOP_CODON:
            result[i] = elem
        else:
            result = result[:i]
            break
    if to_str:
        result = "".join(result)
    return result


def dna_to_protein(dna: str):
    """ Return protein string based on dna string.
        Just a conveyor rna_to_protein(dna_to_rna(dna)), defined for simplifying syntax.
    :param dna: dna string
    :return: protein string
    """
    return rna_to_protein(dna_to_rna(dna))


def gene_to_protein(gene: str, intrones: Union[str, Collection[str]]) -> str:
    """ Return protein for gene with intrones taken into accounts
    :param gene: dna string
    :param intrones: intrones in gene, which will be deleted while generating the protein
    :return: A protein string resulting from transcribing and translating the exons of gene.
    """
    intrones = intrones if not isinstance(intrones, str) else (intrones,)
    for introne in intrones:
        gene = gene.replace(introne, "")
    return dna_to_protein(gene)


def reverse_complement(dna: str):
    """The reverse complement of a DNA string s is the string sc formed by reversing the symbols of s,
    then taking the complement of each symbol (e.g., the reverse complement of "GTCA" is "TGAC").

    Given: A DNA string s of length at most 1000 bp.
    Return: The reverse complement sc of s.
    """
    dna = dna.strip()
    result = [" "] * len(dna)
    for index, letter in enumerate(reversed(dna)):
        result[index] = complement_map[letter]
    return "".join(result)


def fasta_file_iter(path: str):
    """ Return an iterator based on fasta file, so you can use it if for cycle
    Usage:
        profile_matrix = None
        for id, dna in fasta_file_iter("data/test_profile_fasta.txt"):
            print(id, dna)
            profile_matrix = profile(dna, profile_matrix)
            print(profile_matrix)
    :param path: path to file in fasta format
    :return: an iterator to pass to for cycle
    """
    with open(path, "r") as file:
        id = None
        string = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if id:
                    yield id, string
                id = line[1:]
                string = ""
            else:
                string += line
        yield id, string


def fasta_data_do(path, action, result=True):
    """ do some actions for every dna (in fasta format) in file
    DNA strings must be labeled when they are consolidated into a database. A commonly used method of string
        labeling is called FASTA format. In this format, the string is introduced by a line that begins with '>',
        followed by some labeling information. Subsequent lines contain the string itself; the first line to begin
        with '>' indicates the label of the next string. In Rosalind's implementation, a string in FASTA format will
        be labeled by the ID "Rosalind_xxxx", where "xxxx" denotes a four-digit code between 0000 and 9999.
        eg. >Rosalind_6404
            CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
            TCCCACTAATAATTCTGAGG

    :param result: if truely will return result
                   if list then will put all results to it
    """
    result_to_return = None
    for id, dna in fasta_file_iter(path):
        result_to_return = _call_action(action, id, dna, result)
    return result_to_return


def _call_action(action, id, string, result=True):
    if not id:
        return
    if isinstance(result, list):
        result.append(action(id, string))
    elif result:
        result = action(id, string)
    else:
        action(id, string)
    return result


def gc_rate(dna: str, procent=False):
    """ returns rate for G and C in dna
    :param procent:
    """
    c = Counter(dna)
    result = (c["G"] + c["C"]) / len(dna)
    return result * 100 if procent else result


def hamming_distance(p, q):
    """ Compute the Hamming distance between two strings.
    :return: The Hamming distance between these strings.
    >>> hamming_distance("GGGCCGTTGGT", "GGACCGTTGAC")
    3
    """
    result = 0
    for x, y in zip(p, q):
        if x != y:
            result += 1
    return result + abs(len(p) - len(q))


def find_motif(dna:str, motif: str, zero_based=True):
    """ returns indexes of all occurrences of motif in dna.
    :param dna: the string to search in
    :param motif: the substring to search
    :param zero_based: if True will return indexes starting with 1 instead of 0.
    :return: indexes of all occurrences of motif in dna
    """
    def helper_for_non_zero_based(indexes: List[int]):
        if not zero_based:
            return [i + 1 for i in indexes]
        else:
            return indexes

    index = 0
    result = []
    while index >=0:
        index = dna.find(motif, index)
        if index >=0:
            result.append(index)
            index += 1
    return helper_for_non_zero_based(result)


@lru_cache(None)
def rabbits_recurrence(n, k=1):
    if n < 3:
        return 1
    else:
        return rabbits_recurrence(n-1, k) + k * rabbits_recurrence(n-2, k)


@lru_cache(None)
def dying_rabbits(n, months_of_life):
    """
    :return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live
     for m months.
    """
    young = 1
    olds = [0] * (months_of_life - 1)
    for _ in range(n - 1):
        new_young = sum(olds)
        for i in range(1, months_of_life - 1):
            olds[-i] = olds[-i-1]
        olds[0] = young
        young = new_young
    return sum(olds) + young


def calculate_protein_mass(protein: str):
    """ Calculate the standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic
        mass of the corresponding amino acid.
    :param protein: A protein string
    :return: The total weight
    """
    result = 0
    for p in protein:
        result += monoisotopic_mass_table[p]
    return result


def dominant_probability(homozygous_dominant: int, heterozygous: int, homozygous_recessive :int):
    """ Get three positive integers, representing a population containing different genes types,
        return  the probability that two randomly selected mating organisms will produce dominant child.
    :param homozygous_dominant: number of individuals with according genes
    :param heterozygous: number of individuals with according genes
    :param homozygous_recessive: number of individuals with according genes
    :return: the probability that two randomly selected mating organisms will produce an individual possessing a
        dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    """
    d, h, r = homozygous_dominant, heterozygous, homozygous_recessive
    all_ = d + h + r
    result = d * (d + 2 * h + 2 * r - 1) + h * (0.75 * h + r - 0.75)
    result /= all_ * (all_ - 1)
    return result


def profile(dna: Union[list, tuple, str], update: Union[list, None]=None) -> dict:
    """
    Function takes a list of strings DNA as input and returns the profile matrix (as a dictionary of lists).
    :param dna: a list of strings (or just one string) which represent a genome part
    :return: dictionary where keys are A, C, G, T and values are list with their occurrences in patterns
        on that index.
    :example:
    >>> profile(("AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"))
    {'A': [1, 2, 1, 0, 0, 2], 'C': [2, 1, 4, 2, 0, 0], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2]}
    """
    dnas = dna if isinstance(dna, list) or isinstance(dna, tuple) else (dna,)
    k = len(dnas[0])
    update = update if update else {letter: [0] * k for letter in "ACGT"}
    for i in range(k):
        for motif in dnas:
            update[motif[i]][i] += 1
    return update


def consensus(dnas: Union[list, None]=None, precalculated_profile: Union[Dict[str, list], None]=None) -> str:
    """
    Form a consensus string, from the most popular nucleotides in each column of the motif matrix
        (ties are broken arbitrarily). If we select Motifs correctly from the collection of upstream regions,
        then Consensus(Motifs) provides a candidate regulatory motif for these regions.
    :param dnas: A set of kmers.
    :return: A consensus string of dnas.
    :example:
    >>> consensus(("AACGTA","CCCGTT","CACCTT","GGATTA","TTCCGG"))
    'CACCTA'
    """
    k = len(dnas[0]) if dnas else len(precalculated_profile["A"])
    count = precalculated_profile if precalculated_profile else profile(dnas)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def n_reverse_translation(protein: str, modulo: Union[int, None]=None):
    """ returns the total number of different RNA strings from which the protein could have been translated,
     modulo m.
    """
    result = 1
    for p in protein:
        result *= protein_n_codons_table[p]
        if modulo:
            result %= modulo
    return result * protein_n_codons_table[STOP_CODON]


def find_reverse_palindromes(dna: str, min_len: int=4, max_len: int=12, zero_based: bool=True):
    """ A DNA string is a reverse palindrome if it is equal to its reverse complement.
        For instance, GCATGC is a reverse palindrome because its reverse complement is GCATGC.
    :param dna: A DNA string
    :param min_len: minimal length of reversed palindrome to search
    :param max_len: maximal length of reversed palindrome to search
    :param zero_based: if true return indexes starting with 0, or starting with 1, if false
    :return: The position and length of every reverse palindrome in the string having length between min_len and
        max_len.
    """
    def helper_for_non_zero_based(indexes: List[Tuple[int, int]]):
        if not zero_based:
            return [(i + 1, l) for i, l in indexes]
        else:
            return indexes

    length = len(dna)
    result = []
    for i in range(length):
        for l in range(min(min_len, length - i), min(max_len + 1, length - i + 1)):
            if l > max_len or l < min_len:
                continue
            sub_dna = dna[i: i + l]
            if sub_dna == reverse_complement(sub_dna):
                result.append((i, l))
    return helper_for_non_zero_based(result)


def bernul(n, k, p):
    """ returns probability of k occurrences in n Bernoulli trial with probability p
        https://en.wikipedia.org/wiki/Bernoulli_trial
    :param n: number of tests
    :param k: number of successes
    :param p: probability of every success
    :return: probability of k occurrences in n Bernoulli trial
    """
    return comb(n, k) * p ** k * (1 - p) ** (n-k)


def independent_alleles(heterozygous_number: int, generation_number: int) -> float:
    """ http://rosalind.info/problems/lia/
    In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children
        in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism
        having genotype Aa Bb.
    :param heterozygous_number:
    :type generation_number:
    :return: The probability that at least heterozygous_number Aa Bb organisms will belong to the k-th generation of
        Tom's family tree (don't count the Aa Bb mates at each level).
        Assume that Mendel's second law holds for the factors.
    >>> result = independent_alleles(1, 2)
    >>> round(result, 3)
    0.684
    """
    n_child = 2 ** generation_number
    result = 1
    for i in range(0, heterozygous_number):
        result -= bernul(n_child, i, p=1/4)
    return result


def signed_permutation(n: int):
    """ A signed permutation of length n is some ordering of the positive integers {1,2,…,n} in which each integer is
    then provided with either a positive or negative sign (for the sake of simplicity, we omit the positive sign).
    For example, π=(5,−3,−2,1,4) is a signed permutation of length 5.
    :param n: positive integer
    :return: permutations for every digit in [1..n]
    >>> p = signed_permutation(2)
    >>> list(p)
    [[-1, -2], [-1, 2], [1, -2], [1, 2], [-2, -1], [-2, 1], [2, -1], [2, 1]]
    """
    for p in product(permutations(list(range(1, n+1))), product([-1, 1], repeat=n)):
        yield [i * j for i, j in zip(*p)]
