import os
import sys

FLOAT_EQUALITY_ACCURACY = 0.001

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, "test_data_files")))
