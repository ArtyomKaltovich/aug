import random
import string

def random_string(min_len=0, max_len=100):
    return "".join(random.choices(string.ascii_letters, k=random.randint(min_len, max_len)))
