import unittest
import numpy as np
from gillespy2.example_models import Example
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver

class TestAllSolvers(unittest.TestCase):

    solvers = [SSACSolver, BasicODESolver, NumPySSASolver, BasicTauLeapingSolver, BasicTauHybridSolver, CythonSSASolver]
    model = Example()
    results = {}
    labeled_results = {}

    for solver in solvers:
        results[solver] = model.run(solver=solver, show_labels=False, seed=1)
        labeled_results[solver] = model.run(solver=solver, show_labels=True)

    def test_return_type(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.results[solver], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0][0], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0][0][0], np.float))

    def test_return_type_show_labels(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver], list))
            self.assertTrue(isinstance(self.labeled_results[solver][0], dict))
            self.assertTrue(isinstance(self.labeled_results[solver][0]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results[solver][0]['Sp'][0], np.float))

    def test_random_seed(self):
        for solver in self.solvers:
            same_results = self.model.run(solver=solver, show_labels=False, seed=1)
            self.assertTrue(np.array_equal(same_results, self.results[solver]))
            diff_results = self.model.run(solver=solver, show_labels=False, seed=2)
            if solver.name != 'BasicODESolver':
                self.assertFalse(np.array_equal(diff_results, self.results[solver]))
    
    def test_extraneous_args(self):
        for solver in self.solvers:
            print(solver.name)
            with self.assertWarns(Warning):
                model = Example()
                results = model.run(solver=solver, nonsense='ABC')

if __name__ == '__main__':
    unittest.main()
