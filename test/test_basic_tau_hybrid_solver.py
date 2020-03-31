import unittest
import numpy as np
import gillespy2
from example_models import Example
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver
from gillespy2.core.gillespyError import *
from gillespy2.core.results import results

class TestBasicTauHybridSolver(unittest.TestCase):

    model = Example()
    results = model.run(solver=BasicTauHybridSolver, show_labels=False,number_of_trajectories=1)
    labels_results = model.run(solver=BasicTauHybridSolver, show_labels=True,number_of_trajectories=1)

    def test_return_type(self):
        assert(isinstance(self.results, np.ndarray))
        assert(isinstance(self.results[0], np.ndarray))
        assert(isinstance(self.results[0][0], np.ndarray))
        assert(isinstance(self.results[0][0][0], np.float))

    def test_return_type_show_labels(self):
        assert (isinstance(self.labels_results, results))
        assert (isinstance(self.labels_results['Sp'], np.ndarray))
        assert (isinstance(self.labels_results['Sp'][0], np.float))

    def test_add_rate_rule(self):
        species = gillespy2.Species('test_species', initial_value=1, mode='continuous')
        rule = gillespy2.RateRule(species, 'cos(t)')
        self.model.add_species([species])
        self.model.add_rate_rule([rule])
        self.model.run(solver=BasicTauHybridSolver)
        
    def test_add_rate_rule_dict(self):
        species2 = gillespy2.Species('test_species2',initial_value=2, mode ='continuous')
        species3 = gillespy2.Species('test_species3',initial_value=3, mode='continuous')
        rule2 = gillespy2.RateRule(species2, 'cos(t)')
        rule3 = gillespy2.RateRule(species3, 'sin(t)')
        rate_rule_dict = {'rule2' :rule2, 'rule3':rule3}
        self.model.add_species([species2,species3])
        with self.assertRaises(ParameterError):
            self.model.add_rate_rule(rate_rule_dict)
            

if __name__ == '__main__':
    unittest.main()
