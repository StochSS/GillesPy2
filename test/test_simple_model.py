import unittest
import sys, os
sys.path.append(os.path.abspath(os.getcwd()))
import numpy, math, gillespy2
from gillespy2.basic_tau_hybrid_solver import BasicTauHybridSolver
import matplotlib.pyplot as plt


class SimpleHybridModel(gillespy2.Model):
     def __init__(self, parameter_values=None):
            #initialize Model
            gillespy2.Model.__init__(self, name="Simple_Hybrid_Model")
            #Species
            A = gillespy2.Species(name='A', initial_value=0)

            self.add_species([A])

class TestSimpleModel(unittest.TestCase):
    def setUp(self):
        self.model = SimpleHybridModel()
    def test_model_has_species(self):
        species = self.model.get_species('A')
        self.assertEqual(self.model.name, "Simple_Hybrid_Model", msg="{0} is not expected name".format(self.model.name))
        self.assertIsInstance(species, gillespy2.Species, msg='{0} has incorrect type'.format(species))
if __name__ == '__main__':
    unittest.main()