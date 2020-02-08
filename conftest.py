import os
import random
from functools import partial

import numpy
import pytest

FLOAT_EQUALITY_ACCURACY = 0.001
pytest.approx = partial(pytest.approx, rel=FLOAT_EQUALITY_ACCURACY)


@pytest.fixture
def random_seed():
    # TODO: change scope to module and reset it in test dropdown, tearup
    seed = random.randint(0, 10000)
    random.seed(seed)
    numpy.random.seed(seed)
    print("random seed:", seed, end=" ")


@pytest.fixture(scope='module')
def base_data_path():
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    result = os.path.normpath(SCRIPT_DIR)
    result = os.path.join(result, "tests", "test_data_files") + os.sep
    return result
