import random

from aug.seq.seq import edit_distance

if __name__ == '__main__':
    (line1, line2), _ = edit_distance("ab", "ab", reconstract_answer=True)
    print(line1)
    print(line2)
