import urllib.request

from aug.data.fasta import parse_fasta_string


def request_protein_sequence(*ids, source: str = "uniprot"):
    assert source == "uniprot", "currently other sources isn't supported, you can use uniprot only"
    for id in ids:
        content = urllib.request.urlopen(f"https://www.uniprot.org/uniprot/{id}.fasta").read().decode()
        yield from parse_fasta_string(content)
