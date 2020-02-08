from aug.seq.seq import rna_structure_prediction, rna_structure_to_graphviz, dna_to_rna

if __name__ == '__main__':
    rna = "AAACAUGAGGAUUACCCAUGU"
    structure, score = rna_structure_prediction(rna)
    print("score:", score)
    print(rna_structure_to_graphviz(rna, structure))
