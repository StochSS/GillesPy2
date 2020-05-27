import unittest
import numpy as np
import subprocess

from example_models import MichaelisMenten
from gillespy2.core.results import Results, Trajectory
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.variable_ssa_c_solver import VariableSSACSolver
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
from gillespy2.core import gillespyError

class TestPauseResume(unittest.TestCase):
    solvers = [SSACSolver, VariableSSACSolver, BasicODESolver,
               NumPySSASolver, BasicTauLeapingSolver]

    model = MichaelisMenten()
    for sp in model.listOfSpecies.values():
        sp.mode = 'discrete'
    results = {}
    labeled_results = {}
    labeled_results_more_trajectories = {}


    for solver in solvers:
        results[solver] = model.run(solver=solver, show_labels=False)
        labeled_results[solver] = model.run(solver=solver, show_labels=True)

    def test_resume(self):
        model = MichaelisMenten()
        for solver in self.solvers:
            self.results[solver] = model.run(solver=solver, show_labels=False, resume=self.results[solver], t=150)
            self.labeled_results[solver] = model.run(solver=solver, show_labels=True,
                                                     resume=self.labeled_results[solver], t=150)
        for solver in self.solvers:
            self.assertEqual(int(self.results[solver][0][-1][0]),150)
            self.assertEqual(int(self.labeled_results[solver][0]['time'][-1]),150)

    def test_time_fail(self):
        model = MichaelisMenten()
        for solver in self.solvers:
            with self.assertRaises((gillespyError.ExecutionError, gillespyError.SimulationError)):
                self.results[solver] = model.run(solver=solver, show_labels=False, resume=self.results[solver], t=1)
                self.labeled_results = model.run(solver=solver, show_labels=True, resume=self.labeled_results[solver],
                                                 t=1)
