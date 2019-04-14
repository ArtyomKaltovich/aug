import bisect
import random
from copy import copy


def ceil_index(A, key, l=0, r=None):
    r = len(A) if not r else r
    while (r - l > 1):

        m = l + (r - l) // 2
        if (A[m] <= key):
            r = m
        else:
            l = m
    return r


def longest_dec_seq(array):
    # Add boundary case,
    # when array size is one

    if not array:
        return 0

    tail_table = [0] * len(array)
    finish = [-1] * len(array)  # sub seq with len j finished with array[finish[j]]
    previous = [-1] * len(array)

    tail_table[0] = array[0]
    length = 1
    for i in range(1, len(array)):
        if array[i] > tail_table[0]:
            # new smallest value
            tail_table[0] = array[i]
            previous[i] = -1
        elif array[i] <= tail_table[length - 1]:
            # array[i] wants to extend
            # largest subsequence
            tail_table[length] = array[i]
            length += 1
            finish[length] = i
            previous[i] = finish[length - 1]
        else:
            # array[i] wants to be current
            # end candidate of an existing
            # subsequence. It will replace
            # ceil value in tail_table
            pos = ceil_index(tail_table, array[i], r=length)
            tail_table[pos] = array[i]
            previous[i] = finish[pos - 1]

    print(finish)
    print(previous)
    return length, tail_table


def test():
    assert 0 == longest_dec_seq([13, 44, 60, 6, 38, 64, 26, 0, 24, 46, 30, 15, 98])
    assert (5, [(5, 0), (4, 1), (3, 2), (2, 3), (1, 4)]) == longest_dec_seq([5, 4, 3, 2, 1])
    assert (5, [(5, 0), (5, 1), (5, 2), (5, 3), (5, 4)]) == longest_dec_seq([5, 5, 5, 5, 5])
    assert (4, [(5, 0), (4, 2), (4, 3), (2, 4)]) == longest_dec_seq([5, 3, 4, 4, 2])
    assert 0 == longest_dec_seq([])
    assert (1, [(2, 0)]) == longest_dec_seq([2])
    array = [random.randint(0, 100) for _ in range(random.randint(1, 20))]
    print({key + 1: array[key] for key in range(len(array))})
    print(longest_dec_seq(array))


def main():
    n = int(input())
    array = list(map(int, input().split()))
    if len(array) == 10:
        print(*array)
    n, array = longest_dec_seq(array)
    print(n)
    print(*array)


if __name__ == '__main__':
    test()
    #main()