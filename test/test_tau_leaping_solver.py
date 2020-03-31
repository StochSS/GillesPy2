import unittest
import numpy as np
from example_models import Example
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver


class TestBasicTauLeapingSolver(unittest.TestCase):
    model = Example()
    

if __name__ == '__main__':
    unittest.main()
