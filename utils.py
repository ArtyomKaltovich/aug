import itertools
import math
import bisect
from collections import Counter, defaultdict
from functools import lru_cache
from typing import List, Union, Dict, Tuple, Collection, Callable

import numpy as np
from scipy.special import comb

complement_map = {"A": "T", "C": "G", "G": "C", "T": "A"}
START_CODON = "AUG"
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

inverted_monoisotopic_mass = sorted(((value, key) for key, value in monoisotopic_mass_table.items()))


def dna_to_rna(dna: str):
    """An RNA string is a string formed from the alphabet containing 'A', 'C', 'G', and 'U'.

    Given a DNA string t corresponding to a coding strand, its transcribed RNA string u is formed by replacing all
        occurrences of 'T' in t with 'U' in u.
    >>> dna_to_rna("GCAT")
    'GCAU'
    """
    return dna.replace("T", "U")


def start_with_start_codon(rna: str):
    """ helper method to start protein sequence generation from start codon.
    :param rna: rna string
    :return: first position of start codon in rna, or -1 if not presentd
    """
    return rna.find(START_CODON)


def start_with_the_beggining(rna: str):
    """ helper method for processing whole rna sequence starting with 0 index, not START_CODON
    :param rna: rna string
    :return: 0
    """
    return 0


def rna_to_protein(rna: str, to_str=True, start: Union[int, Callable[[str], int]]=start_with_the_beggining):
    """ The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
        Given: An RNA string s corresponding to a strand of mRNA (of length at most 10 kbp).
        Return: The protein string encoded by s.
    :param rna: rna string
    :param to_str: if true result will be returned as string, if false it will be returned as a list
    :param start: start position, can be 0-based int or callable, two variants is already defined:
        * start_with_the_beggining
        * start_with_start_codon
    :return: protein sequence
    >>> seq = "UUUAUGCUUUAA"
    >>> rna_to_protein(seq)
    'FML'
    >>> rna_to_protein(seq, to_str=False)
    ['F', 'M', 'L']
    >>> rna_to_protein(seq, start=start_with_start_codon)
    'ML'
    >>> rna_to_protein(seq, start=6)
    'L'
    """
    pos = start(rna) if callable(start) else start
    if pos < 0:
        return None if not to_str else ""
    else:
        rna = rna[pos:]
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


def all_possible_gene_readings(dna: str):



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
    :param zero_based: if False will return indexes starting with 1 instead of 0.
    :return: indexes of all occurrences of motif in dna
    """

    index = 0
    result = []
    while index >=0:
        index = dna.find(motif, index)
        if index >=0:
            result.append(index)
            index += 1
    return _helper_for_non_zero_based(result, zero_based)


def _helper_for_non_zero_based(indexes: List, zero_based: bool):
    """ Transform indexes based on zero_based.
    :param indexes: list to transform
    :param zero_based: if False indexes will be increased by 1
    :return: list of indexes transformed accordingly to zero_based
    """
    if not zero_based:
        return [i + 1 for i in indexes]
    else:
        return indexes


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
    for p in itertools.product(itertools.permutations(list(range(1, n+1))), itertools.product([-1, 1], repeat=n)):
        yield [i * j for i, j in zip(*p)]


def adjacency_list(fasta_file_path:str, k:int=3, prefixes:Union[str, None]=None,
                   suffixes:Union[str, None]=None) -> Union[List[str], List[Tuple[str, str]]]:
    """ Return adjacency_list of dna in fasta file specified in fasta_file_path.
    For a collection of strings and a positive integer k, return the overlap graph in which each string is represented
        by a node, and string s is connected to string t with a directed edge when there is a length k suffix of s
        that matches a length k prefix of t, as long as s !=t;
        we demand s != t to prevent directed loops in the overlap graph (although directed cycles may be present).
    :param fasta_file_path: path to file with dna in fasta format
    :param k: length of prefixes and suffixes
    :param prefixes: precalculated prefixes or None
    :param suffixes: precalculated suffixes or None
    :return: Overlap graph in form of adjacency list
    """
    prefixes = prefixes or defaultdict(set)
    suffixes = suffixes or defaultdict(set)

    for id, string in fasta_file_iter(fasta_file_path):
        prefixes[string[:k]].add(id)
        suffixes[string[-k:]].add(id)

    result = []
    for suffix, ids in suffixes.items():
        for a in ((start, finish) for start in ids for finish in prefixes[suffix] if start != finish):
            result.append(a)
    return result


def read_fasta(fasta_file_path:str, without_id=False):
    """ read full content of fasta file and return it as a list
    :param fasta_file_path: path to file to read
    :param without_id: if true ids of strings will be omitted
    :return:
    """
    result = []
    for id, string in fasta_file_iter(fasta_file_path):
        if without_id:
            result.append(string)
        else:
            result.append((id, string))
    return result


def dna_probability(dna:str, gc:float, return_log=False) -> float:
    """ For giving dna string and probability of g or c nucleotide return probability of that string or log base 10
    of that probability if return_log is set to True.
    :param dna: dna string
    :param gc: probability of g or c nucleotide (gc-rate)
    :param return_log: set true if you want to get log base 10 of probability
    :return: probability of giving dna string or log base 10 of probability
    >>> result = dna_probability("ACGATACAA", 0.287)
    >>> round(result, 9)
    6.066e-06
    >>> result = dna_probability("ACGATACAA", 0.287, return_log=True)
    >>> round(result, 3)
    -5.217
    """
    at = (1 - gc) / 2.0
    gc /= 2.0

    p = 1
    for l in dna:
        if l in "AT":
            p *= at
        elif l in "CG":
            p *= gc
        else:
            raise ValueError("You should use dna string.")
    if return_log:
        return math.log(p, 10)
    else:
        return p


def find_spliced_motif(dna: str, motif: str, zero_based=True) -> Union[List[int], int]:
    """ Returns the the positions of a subsequence(motif) in the string dna at which the symbols of the subsequence
        appear.
        E.g. the indices of ACG in TATGCTAAGATC can be represented by (2, 5, 9).
    :param dna: dna string
    :param motif: subsequence to search
    :param zero_based: if false will return indexes starting with 1 instead of 0.
    :return: list of indices
    """

    j = 0
    result = []
    for i, l in enumerate(dna):
        if l == motif[j]:
            result.append(i)
            j += 1
        if j >= len(motif):
            break
    else:
        return -1
    return _helper_for_non_zero_based(result, zero_based)


def edit_distance(str1, str2, reconstract_answer=False):
    """ Calculate editing distance between two strings.
    >>> edit_distance("editing", "distance")
    5
    """
    distances = [[0] * (len(str1) + 1) for _ in range(len(str2) + 1)]
    distances[0] = [x for x in range(len(str1) + 1)]
    for x in range(len(str2) + 1):
        distances[x][0] = x
    for i, j in itertools.product(list(range(1, len(str2) + 1)), list(range(1, len(str1) + 1))):
        distances[i][j] = min(distances[i - 1][j] + 1, distances[i][j - 1] + 1,
                              distances[i - 1][j - 1] + (str2[i - 1] != str1[j - 1]))
    if reconstract_answer:
        pass
        # TODO: reconstract answer
    else:
        return distances[-1][-1]


def enumerate_kmers(alphabet: Union[str, List[str]], length: int):
    """ Create generator which will return all words with specified length (k-mers) which can be formed from alphabet.
    :param alphabet:
    :param length: or k in k-mers, length of created words
    :return: all possible words one by one
    >>> result = enumerate_kmers("ab", 3)
    >>> list(result)
    ['aaa', 'aab', 'aba', 'abb', 'baa', 'bab', 'bba', 'bbb']
    >>> result = enumerate_kmers("ABCG", 2)
    >>> list(result)
    ['AA', 'AB', 'AC', 'AG', 'BA', 'BB', 'BC', 'BG', 'CA', 'CB', 'CC', 'CG', 'GA', 'GB', 'GC', 'GG']
    """
    for value in itertools.product(alphabet, repeat=length):
        yield "".join(value)


def string_to_kmers(s: str, k: int) -> List[str]:
    """ Split string S to array of k-mers.
    :param s: string to split
    :param k: length of k-mers
    :return: generator which returns k-mers of split string
    >>> result = string_to_kmers("aaabaa", 2)
    >>> list(result)
    ['aa', 'ab', 'aa']
    >>> result = string_to_kmers("ACGT", 3)
    >>> list(result)
    ['ACG', 'T']
    """
    for i in range(0, len(s), k):
        yield s[i:i + k]


def kmers_composition(dna: str, k: int, alphabet: str = "ACGT"):
    """ For a fixed positive integer k, order all possible k-mers taken from an underlying alphabet lexicographically.
    Then the k-mer composition of a string S can be represented by an array A for which A[m] denotes the number of times
        that the mth k-mer (with respect to the lexicographic order) appears in s.
    :param dna: dna string to represent in k-mer composition
    :param k: length of k-mer
    :param alphabet: alphabet of string
    :return: k-mer composition of dna string
    >>> result = kmers_composition("aaabaa", k=2, alphabet="ab")
    >>> list(result)
    [2, 1, 0, 0]
    """
    dna = Counter(string_to_kmers(dna, k))
    for k_mer in enumerate_kmers(alphabet, k):
        yield dna[k_mer]


def count_kmers(dna: str, k: int, alphabet: str = "ACGT"):
    """ Count number of kmers lexicographically.
    :param dna: dna string to count
    :param k: length of k-mer
    :param alphabet: alphabet of string
    :return: number of k-mer occurs in string in lexicographical order
    """
    c = Counter(dna[i:i + k] for i in range(len(dna) - k + 1))
    result = []
    for k_mer in enumerate_kmers(alphabet, k):
        result.append(c[k_mer])
    return result


def distance_matrix(dnas: Collection[str], metric=hamming_distance, relative=True, as_ndarray=False):
    """ computes matrix distance for string in dnas
    :param dnas: collection of strings
    :param metric: function to calculate distance between two strings
    :param relative: if true distance will be return in 0.0..1.0 interval,
        every item will be divided by size of the biggest string
    :param as_ndarray: if true result will be return as numpy.ndarray
    :return: matrix nxn (where n is length of dnas) where result[i][j] = metric(dnas[i], dnas[j]), possible devided by
        strings size.
    >>> dnas = ["ATTA", "ATTC", "ATTA"]
    >>> distance_matrix(dnas, relative=False)
    [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    """
    n = len(dnas)
    result = [[0] * n for _ in range(n)]
    for pair in itertools.combinations(zip(range(n), dnas), r=2):
        (idx1, dna1), (idx2, dna2) = pair
        distance = metric(dna1, dna2)
        distance = distance / max(len(dna1), len(dna2)) if relative else distance
        result[idx1][idx2] = distance
        result[idx2][idx1] = distance
    if as_ndarray:
        result = np.asarray(result)
    return result


def failure_array(dna: str) -> List[int]:
    """ The failure array of a string is an array P of length n for which P[k] is the length
        of the longest substring s[j:k] that is equal to some prefix s[0:k−j], where j cannot equal 1
        (otherwise, P[k] would always equal k). By convention, P[0]=0.
    :param dna: string to compute failure array from
    :return: computed failure array
    >>> failure_array("CAGCATGGTATCACAGCAGAG")
    [0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 3, 4, 5, 3, 0, 0]
    >>> failure_array("AAAAA")
    [0, 1, 2, 3, 4]
    """
    result = [0] * len(dna)
    for i in range(1, len(dna)):
        for prev in range(result[i - 1], -1, -1):
            if dna[:prev + 1] == dna[i - prev:i + 1]:
                result[i] = prev + 1
                break
    return result


def prefix_spectrum(spectrum: List[float]):
    """ The prefix spectrum of a weighted string is the collection of all its prefix weights.
    :param spectrum: prefix spectrum
    :return: protein sequence, with the same prefix sequence.
        result[i] = aminoacid with the mass equals to spectrum[len(spectrum) - 1] - spectrum[len(spectrum) - 2]
    >>> prefix_spectrum([3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646])
    'WMQS'
    """
    result = []
    it = iter(spectrum)
    next(it)
    m = inverted_monoisotopic_mass
    for pair in zip(spectrum, it):
        mass = pair[1] - pair[0]
        index = bisect.bisect_left(m, (mass, "dummy"))
        index = index - 1 if index == len(m) or\
                             index - 1 > 0 and abs(mass - m[index - 1][0]) < abs(mass - m[index][0]) else index
        result.append(m[index][1])
    return "".join(result)
