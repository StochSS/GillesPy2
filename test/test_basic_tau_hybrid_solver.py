import unittest
import numpy as np
import gillespy2
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver

class TestBasicTauHybridSolver(unittest.TestCase):

    model = Example()
    results = model.run(solver=BasicTauHybridSolver, show_labels=False)
    labels_results = model.run(solver=BasicTauHybridSolver, show_labels=True)

    def test_return_type(self):
        assert(isinstance(self.results, list))
        assert(isinstance(self.results[0], np.ndarray))
        assert(isinstance(self.results[0][0], np.ndarray))
        assert(isinstance(self.results[0][0][0], np.float))

    def test_return_type_show_labels(self):
        assert(isinstance(self.labels_results, list))
        assert(isinstance(self.labels_results[0], dict))
        assert(isinstance(self.labels_results[0]['Sp'], np.ndarray))
        assert(isinstance(self.labels_results[0]['Sp'][0], np.float))

    def test_add_rate_rule(self):
        species = gillespy2.Species('test_species', initial_value=1, mode='continuous')
        rule = gillespy2.RateRule(species, 'cos(t)')
        self.model.add_species([species])
        self.model.add_rate_rule([rule])
        self.model.run(solver=BasicTauHybridSolver)


if __name__ == '__main__':
    unittest.main()
