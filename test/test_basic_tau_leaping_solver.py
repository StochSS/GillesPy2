import unittest
import numpy as np
from example_models import Example
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
from gillespy2.core.results import results

class TestBasicTauLeapingSolver(unittest.TestCase):

    model = Example()
    results = model.run(solver=BasicTauLeapingSolver, show_labels=False,number_of_trajectories=1)
    labels_results = model.run(solver=BasicTauLeapingSolver, show_labels=True,number_of_trajectories=1)

    def test_return_type(self):
        assert(isinstance(self.results, np.ndarray))
        assert(isinstance(self.results[0], np.ndarray))
        assert(isinstance(self.results[0][0], np.ndarray))
        assert(isinstance(self.results[0][0][0], np.float))

    def test_return_type_show_labels(self):
        assert(isinstance(self.labels_results, results))
        assert(isinstance(self.labels_results['Sp'], np.ndarray))
        assert(isinstance(self.labels_results['Sp'][0], np.float))

if __name__ == '__main__':
    unittest.main()
