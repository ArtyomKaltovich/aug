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
        yield from _fasta_structure_iter(file)


def parse_fasta_string(string):
    """ parse one fasta record
    :param string: fasta record in string format (id should start with > and end with \n)
    :return: id and sequence
    >>> parse_fasta_string(">id\\nAUG")
    ('id', 'AUG')
    """
    for r in _fasta_structure_iter(string.splitlines()):
        return r


def _fasta_structure_iter(file):
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