import random
import string


def random_string(min_len=0, max_len=100, alphabet=string.ascii_letters):
    return "".join(random.choices(alphabet, k=random.randint(min_len, max_len)))
