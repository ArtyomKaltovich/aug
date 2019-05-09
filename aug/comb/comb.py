def gen_substrings(alphabet: str, n: int):
    """
    :param alphabet:
    :param n:
    :return:
    >>> list(gen_substrings("DNA", 2))
    ['D', 'DD', 'DN', 'DA', 'N', 'ND', 'NN', 'NA', 'A', 'AD', 'AN', 'AA']
    """
    n = n - 1
    d = [0]
    result = [alphabet[0]]
    while d:
        yield "".join(result)
        if len(d) <= n:
            d.append(0)
            result.append(alphabet[0])
        elif d[-1] < len(alphabet) - 1:
            d[-1] += 1
            result[-1] = alphabet[d[-1]]
        elif d:
            while d and d[-1] >= len(alphabet) - 1:
                d.pop()
                result.pop()
            if d:
                d[-1] += 1
                result[-1] = alphabet[d[-1]]