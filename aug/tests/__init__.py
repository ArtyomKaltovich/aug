import os
import pytest

FLOAT_EQUALITY_ACCURACY = 0.001


@pytest.fixture(scope='module')
def base_data_path():
    SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
    return os.path.normpath(SCRIPT_DIR) + os.sep

