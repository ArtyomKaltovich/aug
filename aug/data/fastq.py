class FastqRecord:
    """ Represents one entity of fastq file

    Parameters
    ----------
    id: str
        id of the record. 1st line (without @ character)
    seq: str
        sequence. 2nd line
    quality_str: str
        quality string
    """

    def __init__(self, id, seq, quality_str):
        self.id = id
        self.seq = seq
        self.quality_str = quality_str

    def __eq__(self, other):
        return other.id == self.id and other.seq == self.seq and other.quality_str == self.quality_str

    def __hash__(self):
        result = hash(self.id)
        for field in ("seq", "quality_str"):
            result += 37 * hash(getattr(self, field))
            result = hash(result)  # to avoid long arithmetic operation
        return result

    def __str__(self):
        return f"FastqRecord(id={self.id}, seq={self.seq}, quality_str={self.quality_str})"

    def __repr__(self):
        return str(self)


def fastq_file_iter(path):
    """ Generator which return

    Parameters
    ----------
    path: str
        Path to fastq file.

    Returns
    -------
    result: generator
        Returns the generator which yields FastqRecord instances
    """
    with open(path) as file:
        id = None
        seq = None
        n = 3
        for i, line in enumerate(file):
            n = i % 4
            if n == 0:
                assert line[0] == "@"
                id = line[1:].strip()
            elif n == 1:
                seq = line.strip()
            elif n == 3:
                yield FastqRecord(id, seq, line.strip())
        assert n == 3, "last record of fastq file seems to be corrupted"
