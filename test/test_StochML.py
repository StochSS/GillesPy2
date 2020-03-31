import unittest
from example_models import Example
from gillespy2.core import Model, StochMLDocument
from gillespy2.core.gillespyError import *
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
import numpy as np


class TestStochML(unittest.TestCase):

    def test_StochML_from_and_to_model(self):
        model = Example()
        stochml = StochMLDocument.from_model(model)
        stochml_model = stochml.to_model('model')
        stochml_model.run(solver=BasicODESolver)
        stochml_model.run(solver=NumPySSASolver)



if __name__ == '__main__':
    unittest.main()
