import itertools
from collections import Counter


class PhenotypeHeredityTable:
    @property
    def table(self):
        if not hasattr(self, "_phenotype_heredity_table") or self._phenotype_heredity_table is None:
            p = ["".join(phenotype) for phenotype in itertools.combinations_with_replacement("Aa", r=2)]
            p = [parents for parents in itertools.combinations_with_replacement(p, r=2)]
            p = {key: Counter("".join(sorted(i)) for i in itertools.product(key[0], key[1]))
                                for key in p}
            # normalize values
            for counter in p.values():
                for key in counter:
                    counter[key] /= 4
            self._phenotype_heredity_table = p
        return self._phenotype_heredity_table

    def dominant_phenotype_prob(self, parents):
        return self.table[parents]["AA"] + self.table[parents]["Aa"]
