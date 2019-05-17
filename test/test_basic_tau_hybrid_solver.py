import unittest
import numpy as np
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver

class TestBasicTauHybridSolver(unittest.TestCase):

    model = Example()
    self.results = model.run(solver=BasicTauHybridSolver, show_labels=False)
    self.labels_results = model.run(solver=BasicTauHybridSolver, show_labels=True)

    def test_return_type(self):
        assert(isinstance(self.results, list))
        assert(isinstance(self.results[0], np.ndarray))
        assert(isinstance(self.results[0]['Sp'], np.ndarray))
        assert(isinstance(self.results[0]['Sp'][0], np.float))

    def test_return_type_show_labels(self):
        assert(isinstance(self.labels_results, list))
        assert(isinstance(self.labels_results[0], dict))
        assert(isinstance(self.labels_results[0]['Sp'], np.ndarray))
        assert(isinstance(self.labels_results[0]['Sp'][0], np.float))


if __name__ == '__main__':
    unittest.main()
