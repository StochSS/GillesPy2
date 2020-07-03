import unittest
from example_models import Example
from gillespy2.core import Model, StochMLDocument
from gillespy2.core.gillespyError import *
from gillespy2 import ODESolver
from gillespy2 import NumPySSASolver


class TestStochML(unittest.TestCase):

    def test_StochML_from_and_to_model(self):
        model = Example()
        stochml = StochMLDocument.from_model(model)
        stochml_model = stochml.to_model('model')
        stochml_model.run(solver=ODESolver)
        stochml_model.run(solver=NumPySSASolver)



if __name__ == '__main__':
    unittest.main()
