import unittest
import numpy as np
from example_models import Example
from gillespy2 import TauLeapingSolver


class TestBasicTauLeapingSolver(unittest.TestCase):
    model = Example()
    

if __name__ == '__main__':
    unittest.main()
