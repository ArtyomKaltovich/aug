import pytest
import numpy
import random

@pytest.fixture
def random_seed():
    seed = random.randint(0, 10000)
    random.seed(seed)
    numpy.random.seed(seed)
    print("random seed:", seed, end=" ")
