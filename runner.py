import random
from functools import partial

from aug.comb.comb import gen_substrings
from aug.data.fasta import *
from aug.data.net import request_protein_sequence, download_from_mg_rast
from aug.seq.seq import *


def read_digits():
    with open("data/rosalind_data.txt", "r") as f:
        line = f.readline().strip()
        n = list(map(int, line.split()))
    return n


def for_fasta():
    result = []
    for id, string in fasta_file_iter(r"C:\Users\artem\AppData\Local\Temp\mgm4443685.3.050.upload.fna"):
        for n in range(30, 60):
            #result = all_possible_gene_transcription(string)
            if string[:n] == string[-n:]:
                result.append(string)
                print(id, string, len(string))
                break
    for r in result:
        print(r, file=open("data/answer.txt", "a"))


def read_file():
    with open("data/rosalind_data.txt", "r") as f:
        for line in f:
            yield line.strip()
        #i, j = read_digits(dna)


def print_statistics(data):
    #print(", ".join("'" + d + "'" for d in data["metagenome_id"]))
    #for id_ in data["metagenome_id"]:
    #    if id_ not in {'mgm4721954.3', 'mgm4667025.3', 'mgm4747906.3', 'mgm4747909.3', 'mgm4443679.3', 'mgm4721953.3', 'mgm4779606.3', 'mgm4667019.3', 'mgm4747911.3', 'mgm4667026.3', 'mgm4624086.3', 'mgm4443683.3', 'mgm4721951.3', 'mgm4721961.3', 'mgm4667023.3', 'mgm4445845.3', 'mgm4721956.3', 'mgm4667035.3', 'mgm4667034.3', 'mgm4572579.3', 'mgm4667020.3', 'mgm4520502.3', 'mgm4721963.3', 'mgm4624083.3', 'mgm4667018.3', 'mgm4667029.3', 'mgm4721958.3', 'mgm4721955.3', 'mgm4667030.3', 'mgm4747910.3', 'mgm4747905.3', 'mgm4594282.3', 'mgm4721952.3', 'mgm4615178.3', 'mgm4667022.3', 'mgm4721962.3', 'mgm4747904.3', 'mgm4594281.3', 'mgm4443682.3', 'mgm4779568.3', 'mgm4667024.3', 'mgm4667028.3', 'mgm4624085.3', 'mgm4721959.3', 'mgm4721957.3', 'mgm4615181.3', 'mgm4667032.3', 'mgm4667027.3', 'mgm4721960.3', 'mgm4624084.3', 'mgm4667033.3', 'mgm4667021.3', 'mgm4615180.3', 'mgm4749478.3', 'mgm4667031.3', 'mgm4601925.3', 'mgm4615179.3', 'mgm4667036.3'}:
    #        print(id_)
    # mgm4443680.3
    # mgm4443681.3
    # mgm4622666.3
    #print(data["biome"].value_counts())
    #print(data["feature"].value_counts())
    #print(data["env_package"].value_counts())
    print(data["material"].value_counts())
    #print(data[data["material"] == "air"][["biome", "feature", "env_package", "material"]])


if __name__ == '__main__':
    #file=open("data/answer.txt", "w")
    #data = read_fasta("data/rosalind_data.txt", without_id=True)
    #print(transition_transversion_ratio(*data))
    #for_fasta()
    #download_from_mg_rast(folder_to_save_=r"C:\Users\artem\Downloads\antarctica_metagenomes",
    #                      #just_description_=True, description_callback_=print_statistics,
    #                      files_to_download_={"cluster.aa90.faa", "cluster.aa90.mapping"},
    #                      path_=["material"],
    #                      country="Antarctica", order="name", limit=300,
    #                      investigation_type="metagenome")
    unite_fasta_in_folder(r"C:\Users\artem\Downloads\antarctica_metagenomes\air", file_name_part_to_id=".")
