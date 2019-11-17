import itertools
import os
from contextlib import contextmanager


def fasta_file_iter(path: str):
    """ Return an iterator based on fasta file, so you can use it in for cycle
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


def fasta(id, sequence):
    """ create fasta record.
    :param id: str
    :param sequence: str
    :return: string in fasta format:
        >id\n
        seq\n
    """
    return f">{id}\n{sequence}\n"


def unite_fasta(union_path, path, *paths, file_name_part_to_id=None, id_sep=".", append_sep="__"):
    """ unite two fasta files
    Read fasta files and concat the accordingly:
        * remove dublicates (records with the same id and sequence)
        * change id if record with the same id, but different sequence is already presented
    :param union_path:
    :param path: path to the first fast file
    :param paths: path to other fasta files, can be omitted, in such case resulted file will be the same as in path,
        but with processed duplications and optionally with file_id add.
    :param file_name_part_to_id: str, if presented the file name before this characters will be appended to fasta ids.
        These appends don't affect duplications search
    :param append_sep: str, string to add to current id if current is duplicated. After append_sep index will be added.
    :return:
    """

    with open(union_path, "w") as union_file:
        records = {}
        for fasta_file in itertools.chain([path], paths):
            file_id = _get_file_id(fasta_file, file_name_part_to_id, id_sep) if file_name_part_to_id else ""
            for id, seq in fasta_file_iter(fasta_file):
                presented = records.get(id, None)
                if presented is not None:
                    if presented != seq:
                        id = _gen_new_id(id, records, seq, append_sep)
                    if presented == seq or id is None:  # while generating new id, duplicate can be found
                        continue
                union_file.write(fasta(file_id + id, seq))
                records[id] = seq


def _gen_new_id(id: str, records, seq, append_sep):
    before, sep, index = id.rpartition(append_sep)
    while True:
        if index and index.isdecimal():
            i = int(index) + 1
            result = before + append_sep + str(i)
        else:
            result = id + append_sep + "1"
        existed = records.get(result, None)
        if existed is None:
            return result
        elif existed == seq:
            return  # return None if already presented
        before, sep, index = result.rpartition(append_sep)


def _get_file_id(path: str, file_name_part_to_id: str, id_sep: str):
    folder, sep, file = path.rpartition(os.sep)
    file_id, _, _ = file.partition(file_name_part_to_id)
    return f"{file_id}{id_sep}"


@contextmanager
def _with_dir(folder: str):
    prev = os.getcwd()
    os.chdir(folder)
    yield
    os.chdir(prev)


def is_fasta_extension(path, fasta_extensions={"fasta", "faa"}):
    _, _, ext = path.rpartition(".")
    return ext in fasta_extensions


def unite_fasta_in_folder(folder: str, fasta_extensions={"fasta", "faa"}, file_name_part_to_id=None, id_sep=".",
                          append_sep="__"):
    folder = folder if not folder.endswith(os.sep) else folder[:]
    for dirpath, dirnames, filenames in os.walk(folder):
        filenames = [name for name in filenames if is_fasta_extension(name, fasta_extensions)]
        with _with_dir(dirpath):
            print(dirpath + ".fasta", *filenames, file_name_part_to_id, id_sep, append_sep)
            unite_fasta(dirpath + ".fasta", *filenames, file_name_part_to_id=file_name_part_to_id, id_sep=id_sep,
                        append_sep=append_sep)
