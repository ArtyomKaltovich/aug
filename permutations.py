import itertools
import math


if __name__ == '__main__':
    n = 5
    print(math.factorial(n))
    for i in itertools.permutations([i for i in range(1, n + 1)]):
        print(*i)
