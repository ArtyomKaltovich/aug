class Sets:
    """ The system of sets. https://en.wikipedia.org/wiki/Disjoint-set_data_structure
    Implement as an array size n, where:
        0..i..n - elements of sets,
        array[i] - id of set contains element i
        if array[i] == i then this element is a root
    Supported operations:
        find - return set (the biggest one), current element belongs to
        unite(a, b) - unite two sets: a and b
    Fields:
        parents: list - array with ids of parents
        n_disjoint: int - number of disjoint sets in system
    >>> sets = Sets(3)
    >>> sets.unite(0, 1)
    >>> sets.unite(0, 2)
    >>> print(sets.n_disjoint)
    1
    >>> sets = Sets(4)
    >>> sets.unite(0, 1)
    >>> sets.unite(3, 2)
    >>> print(sets.n_disjoint)
    2
    """
    def __init__(self, n):
        self.parents = [i for i in range(n)]
        self.n_disjoint = n

    def find(self, a):
        parent = self.parents[a]
        while parent != a:
            old = parent
            parent = self.parents[parent]
            self.parents[a] = parent
            a = old
        return parent

    def unite(self, a, b):
        parent_a = self.find(a)
        parent_b = self.find(b)
        if parent_a != parent_b:
            self.n_disjoint -= 1
            self.parents[parent_b] = parent_a
