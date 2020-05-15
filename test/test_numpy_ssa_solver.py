import unittest
import numpy as np
from example_models import Example
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver


class TestNumPySSASolver(unittest.TestCase):
    model = Example()
    

if __name__ == '__main__':
    unittest.main()
