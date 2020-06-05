import unittest
import numpy as np
import subprocess

from example_models import MichaelisMenten, Oregonator
from gillespy2.core.results import Results, Trajectory
from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.variable_ssa_c_solver import VariableSSACSolver
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
from gillespy2.core import gillespyError
import subprocess
import signal
import time
import os


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

    def test_pause(self):
        args = [['python3', 'pause_model.py', 'NumPySSASolver'], ['python3', 'pause_model.py', 'BasicTauLeapingSolver'],
                ['python3', 'pause_model.py', 'BasicODESolver']]
        for arg in args:
            p = subprocess.Popen(arg, preexec_fn=os.setsid, stdout=subprocess.PIPE)
            time.sleep(2)
            os.kill(p.pid, signal.SIGINT)
            out, err = p.communicate()
            #End time for Oregonator is 5. If indexing into a numpy array using the form:
            #results[0][-1][0] (where .run(show_labels=False), this index being the last time in the index
            #One would get an output of "5.0" before converting it to an int. Hence, assert time != 5.0 rather than 5.
            self.assertFalse(out.decode('utf-8').rstrip() == '5.0')

        solvers = [VariableSSACSolver,SSACSolver]
        #For the C solvers, timeouts behave identical to a keyboard interrupt, and would return the same data, if
        #one was to KeyBoardInterrupt, at the same time a timeout ended. This is because timeouts in C solvers
        #and KeyBoardInterrupt send the same signal to the subprocess, the only difference is KeyBoardInterrupt is
        #manual, whereas timeout is a set variable
        for solver in solvers:
            model = Oregonator()
            results = model.run(solver=solver,timeout=1,show_labels=False)
            self.assertFalse(results[0][-1][0] == '5.0')



