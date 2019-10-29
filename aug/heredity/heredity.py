import copy
import heapq

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
            pass

        def recalculate_tree_score(self, dist_matrix, elem, left, right):
            left_size = 1 if type(left) is str else left.size
            right_size = 1 if type(right) is str else right.size
            return (dist_matrix[(str(elem), str(left))] * left_size + dist_matrix[(str(elem), str(right))] * right_size)\
                   / (left_size + right_size)

    def weighted(self):
        return PhylogenyTree._WeightedElem

    def unweighted(self):
        return PhylogenyTree._UnweightedElem

    def __init__(self, dist_matrix, method=weighted, sep="-"):
        """
        :param dist_list:
        :param name_to_index:
        :param method:
        :return:

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
        unprocessed = set()
        for l, r in dist_matrix:
            unprocessed.add(l)
            unprocessed.add(r)
        heapq.heapify(dist_list)
        dist_matrix.update({(r, l): score for (l, r), score in dist_matrix.items()})
        return dist_matrix, dist_list, unprocessed

    def __str__(self):
        return self._tree.to_newick()
