import unittest
import numpy as np

import gillespy2
from example_models import Example, Oregonator
from gillespy2.core.results import EnsembleResults, Results
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver


class TestAllSolvers(unittest.TestCase):

    solvers = [SSACSolver, BasicODESolver, NumPySSASolver, BasicTauLeapingSolver, BasicTauHybridSolver]

    model = Example()
    for sp in model.listOfSpecies.values():
        sp.mode = 'discrete'
    results = {}
    labeled_results = {}
    labeled_results_more_trajectories = {}

    for solver in solvers:
        results[solver] = model.run(solver=solver, show_labels=False, seed=1)
        labeled_results[solver] = model.run(solver=solver, show_labels=True,number_of_trajectories=1)
        labeled_results_more_trajectories[solver] = model.run(solver=solver, show_labels=True,number_of_trajectories=2)

    def test_instantiated(self):
        for solver in self.solvers:
            self.model.run(solver=solver())

    def test_return_type(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.results[solver], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0][0], np.ndarray))
            self.assertTrue(isinstance(self.results[solver][0][0][0], np.float))

    def test_return_type_show_labels(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver], Results))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'][0], np.float))

            with self.assertWarns(Warning):
                self.assertTrue(isinstance(self.labeled_results[solver][0], Results))

            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver], EnsembleResults))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0], Results))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'][0], np.float))


    def test_random_seed(self):
        for solver in self.solvers:
            same_results = self.model.run(solver=solver, show_labels=False, seed=1)
            self.assertTrue(np.array_equal(same_results, self.results[solver]))
            if solver.name == 'BasicODESolver': continue
            diff_results = self.model.run(solver=solver, show_labels=False, seed=2)
            self.assertFalse(np.array_equal(diff_results, self.results[solver]))
    
    def test_extraneous_args(self):
        for solver in self.solvers:
            with self.assertLogs(level='WARN'):
                model = Example()
                results = model.run(solver=solver, nonsense='ABC')

    def test_timeout(self):
        for solver in self.solvers:
            with self.assertLogs(level='WARN'):
                model = Oregonator()
                model.timespan(np.linspace(0, 1000000, 101))
                results = model.run(solver=solver, timeout=1)

if __name__ == '__main__':
    unittest.main()
