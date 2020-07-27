import unittest
import numpy as np

import gillespy2
from example_models import Example, Oregonator, MichaelisMenten
from gillespy2.core.results import Results, Trajectory
from gillespy2 import SSACSolver
from gillespy2 import VariableSSACSolver
from gillespy2 import ODESolver
from gillespy2 import NumPySSASolver
from gillespy2 import TauLeapingSolver
from gillespy2 import TauHybridSolver


class TestAllSolvers(unittest.TestCase):

    solvers = [SSACSolver, VariableSSACSolver, ODESolver, NumPySSASolver, TauLeapingSolver, TauHybridSolver]

    model = Example()
    for sp in model.listOfSpecies.values():
        sp.mode = 'discrete'
    results = {}
    labeled_results = {}
    labeled_results_more_trajectories = {}

    for solver in solvers:
        labeled_results[solver] = model.run(solver=solver, number_of_trajectories=1,seed=1)
        labeled_results_more_trajectories[solver] = model.run(solver=solver, number_of_trajectories=2)

    def test_instantiated(self):
        for solver in self.solvers:
            self.model.run(solver=solver())

    def test_to_array(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver].to_array()[0], np.ndarray))

    def test_return_type_show_labels(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver], Results))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'][0], np.float))

            self.assertTrue(isinstance(self.labeled_results[solver][0], Trajectory))

            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver], Results))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0], Trajectory))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'][0], np.float))


    def test_random_seed(self):
        for solver in self.solvers:
            same_results = self.model.run(solver=solver, seed=1)
            compare_results = self.model.run(solver=solver,seed=1)
            self.assertTrue(np.array_equal(same_results.to_array(), compare_results.to_array()))
            if solver.name == 'ODESolver': continue
            diff_results = self.model.run(solver=solver, seed=2)
            self.assertFalse(np.array_equal(diff_results.to_array(),same_results.to_array()))
    
    def test_random_seed_unnamed_reactions(self):
        model = self.model
        k2 = gillespy2.Parameter(name='k2', expression=1.0)
        model.add_parameter([k2])
        unnamed_rxn = gillespy2.Reaction(reactants={}, products={'Sp':1}, rate=k2)
        model.add_reaction(unnamed_rxn)
        for solver in self.solvers:
            same_results = self.model.run(solver=solver, seed=1)
            compare_results = self.model.run(solver=solver,seed=1)
            self.assertTrue(np.array_equal(same_results.to_array(), compare_results.to_array()))
            if solver.name == 'ODESolver': continue
            diff_results = self.model.run(solver=solver, seed=2)
            self.assertFalse(np.array_equal(diff_results.to_array(),same_results.to_array()))

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

    def test_basic_solver_import(self):
        from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
        from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
        from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver
        model = MichaelisMenten()
        results1 = model.run(solver=BasicTauLeapingSolver)
        self.assertTrue(results1[0].solver_name == 'TauLeapingSolver')

        results2 = model.run(solver=BasicODESolver)
        self.assertTrue(results2[0].solver_name == 'ODESolver')

        results3 = model.run(solver=BasicTauHybridSolver)
        self.assertTrue(results3[0].solver_name == 'TauHybridSolver')


if __name__ == '__main__':
    unittest.main()
