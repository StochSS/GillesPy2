import unittest
import numpy as np
from gillespy2.example_models import Example
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver


class TestBasicTauLeapingSolver(unittest.TestCase):

    model = Example()
    results = model.run(solver=BasicTauLeapingSolver, show_labels=False)
    labels_results = model.run(solver=BasicTauLeapingSolver, show_labels=True)

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

    # Issue #160
    def test_kwargs_throw_error(self):
        solver = BasicTauLeapingSolver

        with self.assertRaises(TypeError):
            solver.run(model=Example(), max_steps=0)

if __name__ == '__main__':
    unittest.main()
