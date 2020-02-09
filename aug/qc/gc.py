import seaborn as sns

from functools import singledispatch
from typing import List

from aug.data.fasta import fasta_file_iter
from aug.seq.seq import gc_rate


@singledispatch
def gc_percent_hist(data):
    """ Calculate the gc rate (in percent) for every record

    Parameters
    ----------
    data: Union[str, List[str]]
        Either path to fasta file or a list of sequences

    Returns
    -------
    result: List[int]
        List of ints, where result[i] = number of record with gc rate i%

    Examples
    --------

    >>> gc_percent_hist("sss")
    'sss'
    >>> gc_percent_hist(["AGC", "AAA", "CTA"])
    [1, ..., 1, ..., 1, ...]
    """
    result = [0] * 101  # 0-100
    for s in data:
        percent = _gc_percent(s)
        result[percent] += 1
    return result


@gc_percent_hist.register(str)
def gc_percent_file(data):
    result = [0] * 101  # 0-100
    for id, seq in fasta_file_iter(data):
        result[_gc_percent(seq)] += 1
    return result


def _gc_percent(s):
    return int(round(gc_rate(s) * 100))


@singledispatch
def gc_plot(data, **kwargs):
    """ Draw gc content histogram using seaborn library.
    Plotting on matplotlib plt object. You can then directly show `savefig`
    or `show` method of it to save of show plot.

    Parameters
    ----------
    data
    kwargs

    Returns
    -------
    None:
    """
    hist = gc_percent_hist(data)
    sns.lineplot(x=range(0, 101), y=hist, **kwargs)
