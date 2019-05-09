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