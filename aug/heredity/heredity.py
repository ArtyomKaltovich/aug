import copy
import heapq
import itertools
from typing import Dict, Tuple

import numpy as np
from scipy.stats import chi2

from aug.heredity.Phenotype import PhenotypeHeredityTable


def n_expected_dominant_phenotype(n_parents, n_children=1):
    """
    :param n_parents:
    :param n_children:
    :return:
    >>> n_expected_dominant_phenotype([1, 0, 0, 1, 0, 1], n_children=2)
    3.5
    """
    assert 6 == len(n_parents)
    result = 0
    phenotype = PhenotypeHeredityTable()
    for n, pair in zip(n_parents, (("AA", "AA"), ("AA", "Aa"), ("AA", "aa"), ("Aa", "Aa"), ("Aa", "aa"), ("aa", "aa"))):
        result += n * phenotype.dominant_phenotype_prob(pair) * n_children
    return result


def test_hwe(n_homo_dominant, n_hetero, n_homo_recessive) -> float:
    """ Statistical test for Hardy–Weinberg equilibrium

    Parameters
    ----------
    n_homo_dominant: float, int
        Number of dominant homozygous (usually labeled as AA)
    n_hetero: float, int
        Number of heterozygous (Aa)
    n_homo_recessive: float, int
        Number of recessive homozygous (Aa)

    Returns
    -------
    result: float
        Probability of such distribution if Hardy–Weinberg equilibrium is true

    Examples
    --------
    >>> test_hwe(313, 102, 85)
    1.0
    """
    actual = [n_homo_dominant, n_hetero, n_homo_recessive]
    n = sum(actual)
    p = (n_homo_dominant * 2 + n_hetero) / n / 2
    q = (n_homo_recessive * 2 + n_hetero) / n / 2
    expected = np.array([p ** 2, 2 * p * q, q ** 2])
    expected *= n
    s = (actual - expected) ** 2 / expected
    s = sum(s)
    return 1 - chi2(2).cdf(s)


class PhylogenyTree:
    class _WeightedElem:
        def __init__(self, left=None, right=None, score=None, sep="-"):
            self.left = left
            self.right = right
            self.score = score
            self.name = f"{self.left}{sep}{self.right}"

        def to_newick(self):
            left, left_score = (self.left, self.score / 2) if type(self.left) is str else\
                (self.left.to_newick(), (self.score - self.left.score) / 2)
            right, right_score = (self.right, self.score / 2) if type(self.right) is str else\
                (self.right.to_newick(), (self.score - self.right.score) / 2)
            return f"({left}:{left_score}, {right}:{right_score})"

        def recalculate_tree_score(self, dist_matrix, elem, left, right):
            return (dist_matrix[(str(elem), str(left))] + dist_matrix[(str(elem), str(right))]) / 2

        def __str__(self):
            return self.name

        def __lt__(self, other):
            # for heapq
            return self.name.__lt__(other)

        def __gt__(self, other):
            # for heapq
            return self.name.__gt__(other)

    class _UnweightedElem(_WeightedElem):
        def __init__(self, left=None, right=None, score=None, sep="-"):
            super().__init__(left, right, score, sep)
            self.size = left.size if isinstance(left, PhylogenyTree._UnweightedElem) else 1
            self.size += right.size if isinstance(right, PhylogenyTree._UnweightedElem) else 1

        def recalculate_tree_score(self, dist_matrix, elem, left, right):
            left_size = 1 if type(left) is str else left.size
            right_size = 1 if type(right) is str else right.size
            return (dist_matrix[(str(elem), str(left))] * left_size + dist_matrix[(str(elem), str(right))] * right_size)\
                   / (left_size + right_size)

    class _NeighbourJoining(_WeightedElem):
        def __init__(self, left=None, right=None, score=None, sep="-"):
            super().__init__(left, right, score, sep)

        def to_newick(self):
            left = self.left if type(self.left) is str else self.left.to_newick()
            right = self.right if type(self.right) is str else self.right.to_newick()
            return f"({left}:{self.score[0]}, {right}:{self.score[1]})"

    def weighted(self):
        return PhylogenyTree._WeightedElem

    def unweighted(self):
        return PhylogenyTree._UnweightedElem

    def neighbour_joining(self):
        return PhylogenyTree._NeighbourJoining

    def __init__(self, dist_matrix: Dict[Tuple[str, str], float], method=weighted, sep="-"):
        """ Builds Phylogeny Tree
        Can be translated to newick format by str function
        :param dist_matrix: matrix of distances
        :param method: method to build tree
        :param sep: separator between names in united tree

        >>> dist_matrix = {("A", "B"): 16, ("A", "C"): 16, ("A", "D"): 10,
        ...                ("B", "C"): 8, ("B", "D"): 8,
        ...                ("C", "D"): 4}
        >>> tree = PhylogenyTree(dist_matrix)
        >>> str(tree)
        '(A:7.25, (B:4.0, (C:2.0, D:2.0):2.0):3.25)'
        >>> dist_matrix = {("A", "B"): 5, ("A", "C"): 4, ("A", "D"): 7, ("A", "E"): 6, ("A", "F"): 8,
        ...                ("B", "C"): 7, ("B", "D"): 10, ("B", "E"): 9, ("B", "F"): 11,
        ...                ("C", "D"): 7, ("C", "E"): 6, ("C", "F"): 8,
        ...                ("D", "E"): 5, ("D", "F"): 9,
        ...                ("E", "F"): 8}
        >>> tree = PhylogenyTree(dist_matrix)
        >>> str(tree)
        '(F:4.5, ((D:2.5, E:2.5):1.5, (B:3.0, (A:2.0, C:2.0):1.0):1.0):0.5)'

        >>> dist_matrix = {("A", "B"): 16, ("A", "C"): 16, ("A", "D"): 10,
        ...                ("B", "C"): 8, ("B", "D"): 8,
        ...                ("C", "D"): 4}
        >>> tree = PhylogenyTree(dist_matrix, method=PhylogenyTree.unweighted)
        >>> str(tree)
        '(A:7.0, (B:4.0, (C:2.0, D:2.0):2.0):3.0)'
        >>> dist_matrix = {("A", "B"): 5, ("A", "C"): 4, ("A", "D"): 7, ("A", "E"): 6, ("A", "F"): 8,
        ...                ("B", "C"): 7, ("B", "D"): 10, ("B", "E"): 9, ("B", "F"): 11,
        ...                ("C", "D"): 7, ("C", "E"): 6, ("C", "F"): 8,
        ...                ("D", "E"): 5, ("D", "F"): 9,
        ...                ("E", "F"): 8}
        >>> tree = PhylogenyTree(dist_matrix, method=PhylogenyTree.unweighted)
        >>> str(tree)
        '(F:4.4, ((D:2.5, E:2.5):1.25, (B:3.0, (A:2.0, C:2.0):1.0):0.75):0.6500000000000004)'
        """
        if method == PhylogenyTree.weighted or method == PhylogenyTree.unweighted:
            # using heapq to decrease complexity of min search
            # 1. save all values to heapq
            # 2. select min
            # 3. add dists for united values to heapq
            # 4. until there are elements in heapq go to step 2
            # On step 2 there can be invalid value selected (distance to element, that was already united)
            # the algorithm will keep all united values in set, and if selected minimum is for invalid elements,
            # it will be through away.
            # So with dist matrix we should find min O(N^2) and update row\column O(N) on every step.
            # Now the algorithm will discard invalid values (on each step the can be only O(N) values).
            # While inserting is O(logN^2) and we will insert N value, so complexity of one step will be O(NlogN^2)
            dist_matrix, dist_list, unprocessed = self._get_initial_values_for_tree(dist_matrix)
            method = method(self)
            self._tree = method()
            while len(unprocessed) > 1:
                left, right, score = self._get_next_tree_elem(dist_list, unprocessed)
                self._tree = method(left, right, score, sep)
                unprocessed.discard(left)
                unprocessed.discard(right)
                new_elem = self._tree
                for elem in unprocessed:
                    self._update_tree_score(dist_matrix, dist_list, elem, new_elem, left, right)
                unprocessed.add(new_elem)
        elif method == PhylogenyTree.neighbour_joining:
            unprocessed = self._get_names_set_from_dist_matrix(dist_matrix)
            name_to_index = {}
            index_to_name = []
            for i, key in enumerate(unprocessed):
                name_to_index[key] = i
                index_to_name.append(key)
            unprocessed = set(range(len(index_to_name)))
            matrix = [[0] * len(name_to_index) for _ in range(len(name_to_index))]
            for (left, right), score in dist_matrix.items():
                i, j = name_to_index[left], name_to_index[right]
                matrix[i][j] = score
                matrix[j][i] = score

            method = method(self)
            self._tree = method()

            for n in range(len(matrix), 1, -1):
                dist_sums = [sum(row) for row in matrix]
                min_ = [None, None, None, None, None]  # score, i, j, left_score, right_score
                for i, j in itertools.combinations(unprocessed, 2):
                    score = matrix[i][j] - (dist_sums[i] + dist_sums[j] - 2 * matrix[i][j]) / (n - 2) if n > 2 else matrix[i][j]
                    if min_[0] is None or score < min_[0]:
                        min_[0] = score
                        min_[1] = i
                        min_[2] = j
                        temp = (dist_sums[i] - dist_sums[j]) / (n - 2) if n > 2 else 0
                        min_[3] = (matrix[i][j] + temp) / 2
                        min_[4] = (matrix[i][j] - temp) / 2
                left = min_[1]
                right = min_[2]
                self._tree = method(index_to_name[left], index_to_name[right], (min_[3], min_[4]), sep)

                name_to_index[str(self._tree)] = left
                index_to_name[left] = self._tree
                unprocessed.discard(right)
                m = copy.deepcopy(matrix)
                for k in unprocessed - {left}:
                    score = (matrix[left][k] + matrix[right][k] - matrix[left][right]) / 2
                    m[left][k] = score
                    m[k][left] = score
                matrix = m
                for k in range(len(matrix)):
                    matrix[right][k] = 0
                    matrix[k][right] = 0
                pass
        else:
            raise ValueError("Unsupported method!")

    def _update_tree_score(self, dist_matrix, dist_list, elem, new_elem, left, right):
        elem_str = str(elem)
        score = self._tree.recalculate_tree_score(dist_matrix, elem, left, right)
        heapq.heappush(dist_list, (score, (elem, new_elem)))
        new_elem_str = str(new_elem)
        dist_matrix[(elem_str, new_elem_str)] = score
        dist_matrix[(new_elem_str, elem_str)] = score

    def _get_next_tree_elem(self, dist_list, unprocessed):
        left, right, score = None, None, None
        while not left or left not in unprocessed or right not in unprocessed:
            score, (left, right) = heapq.heappop(dist_list)
        return left, right, score

    def _get_initial_values_for_tree(self, dist_matrix):
        dist_matrix = copy.copy(dist_matrix)
        dist_list = list((value, key) for key, value in dist_matrix.items())
        unprocessed = self._get_names_set_from_dist_matrix(dist_matrix)
        heapq.heapify(dist_list)
        dist_matrix.update({(r, l): score for (l, r), score in dist_matrix.items()})
        return dist_matrix, dist_list, unprocessed

    def _get_names_set_from_dist_matrix(self, dist_matrix):
        unprocessed = set()
        for l, r in dist_matrix:
            unprocessed.add(l)
            unprocessed.add(r)
        return unprocessed

    def __str__(self):
        return self._tree.to_newick()
