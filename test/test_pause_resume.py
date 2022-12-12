# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import numpy as np
import signal
import threading
import multiprocessing
from example_models import MichaelisMenten, Oregonator
from gillespy2.core.model import Model
from gillespy2.core.results import Results, Trajectory
from gillespy2.core import Species
from gillespy2 import NumPySSASolver
from gillespy2 import ODESolver
from gillespy2 import TauLeapingSolver
from gillespy2 import TauHybridSolver
from gillespy2 import SSACSolver
from gillespy2 import ODECSolver
from gillespy2 import TauLeapingCSolver
from gillespy2 import TauHybridCSolver
from gillespy2.core import gillespyError
import os


def timeout_process(t, pid):
    """
    Delegate function executed by the `multiprocessing` package,
    which sends a keyboard interrupt to the specified process.
    The test runner is expected to execute independently,
    reacting to the interrupt sent by this background process.
    """
    os_signal = signal.SIGINT if os.name != "nt" else signal.CTRL_C_EVENT
    def timeout_kill():
        os.kill(pid, os_signal)

    while True:
        timer = threading.Timer(t, timeout_kill)
        try:
            timer.start()
            timer.join()
        except:
            timer.cancel()

def hybrid_model(model: Model, mode: str = 'discrete'):
    """
    Accepts and returns a model after modifying all species to use the specified mode.
    Used as a clean way to declare the same model for multiple hybrid solver configs.
    """
    for spec in model.listOfSpecies.values():
        spec.mode = mode
    return model

def model_time(model: Model, time: float, increment: float = 1.0):
    """
    Accepts and returns a model after modifying its timespan to use the given time range.
    Used as a clean way to declare the same model with different timespans.
    """
    model.timespan(np.linspace(0, time, int(1 + (time / increment))))
    return model


class TestPauseResume(unittest.TestCase):
    py_solvers = [
        ODESolver,
        NumPySSASolver,
        TauLeapingSolver,
        # TauHybridSolver,
    ]
    c_solvers = [
        SSACSolver,
        ODECSolver,
        TauLeapingCSolver,
        TauHybridCSolver,
    ]
    solvers = [
        *py_solvers,
        *c_solvers,
    ]

    # Models used for testing different solvers' interrupts.
    interrupt_models = {
        "ode": model_time(Oregonator(), 10_000, 0.01),
        "ssa": Oregonator(),
        "discrete": hybrid_model(Oregonator(), 'discrete'),
        "continuous": hybrid_model(Oregonator(), 'continuous'),
    }

    # Match each solver class to a preferred model to construct it with.
    # References the keys in the `interrupt_models` table.
    interrupt_solvers = {
        ODESolver: ["ode"],
        NumPySSASolver: ["ssa"],
        TauLeapingSolver: ["ssa"],
        ODECSolver: ["ode"],
        SSACSolver: ["ssa"],
        TauLeapingCSolver: ["ssa"],
        # TauHybridSolver: ["continuous", "discrete"],
        TauHybridCSolver: ["continuous", "discrete"],
    }

    model = MichaelisMenten()
    labeled_results = {}
    pre_built_solvers = []
    labeled_results_more_trajectories = {}

    @staticmethod
    def setUpClass():
        for sp in TestPauseResume.model.listOfSpecies.values():
            sp.mode = 'discrete'
        for solver_class in TestPauseResume.solvers:
            solver = solver_class(model=TestPauseResume.model)
            TestPauseResume.labeled_results[solver.name] = solver.run()

        # Construct pre-built solvers for interrupt tests
        # These are intentionally designed to run ~30s and are hand-picked for each solver.
        for solver_class, model_type_list in TestPauseResume.interrupt_solvers.items():
            for model_type in model_type_list:
                model = TestPauseResume.interrupt_models.get(model_type)
                solver = solver_class(model=model)
                TestPauseResume.pre_built_solvers.append(solver)

    def test_altered_model_failure(self):
        model = MichaelisMenten()
        for solver_class in self.solvers:
            with self.subTest(solver=solver_class.name):
                solver = solver_class(model=model)
                tmpResults = model.run(solver=solver)
                with self.assertRaises(gillespyError.SimulationError):
                    sp1 = Species('sp2',initial_value=5)
                    model.add_species(sp1)
                    solver = solver_class(model=model)
                    tmpResults = model.run(solver=solver,resume=tmpResults,t=150)
                model.delete_species('sp2')

    def test_resume(self):
        model = self.model
        for solver in self.solvers:
            solver = solver(model=model)
            self.labeled_results[solver.name] = model.run(solver=solver, show_labels=True,
                                                    resume=self.labeled_results[solver.name], t=150)
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                self.assertEqual(int(self.labeled_results[solver.name][0]['time'][-1]),150)

    def test_time_fail(self):
        model = self.model
        for solver in self.solvers:
            with self.subTest(solver=solver.name), self.assertRaises((gillespyError.ExecutionError, gillespyError.SimulationError)):
                solver = solver(model=model)
                self.labeled_results[solver.name] = model.run(solver=solver, show_labels=True, resume=self.labeled_results[solver.name], t=1)

    def test_pause(self):
        for solver in self.pre_built_solvers:
            print("Testing Solver: ", solver.name)
            with self.subTest("Send SIGINT during C++ solver execution", solver=solver.name):
                try:
                    repeater = multiprocessing.Process(target=timeout_process, args=(1.0, os.getpid()))
                    repeater.start()
                    result = solver.run()
                except KeyboardInterrupt:
                    self.fail("solver.run() did not properly handle KeyboardInterrupt exception")
                finally:
                    repeater.terminate()
                max_time = result["time"].max()
                self.assertGreater(max_time, 0.0)
                self.assertLess(max_time, solver.model.tspan[-1])

    def test_timeout(self):
        # For the C solvers, timeouts behave identical to a keyboard interrupt, and would return the same data, if
        # one was to KeyBoardInterrupt, at the same time a timeout ended. This is because timeouts in C solvers
        # and KeyBoardInterrupt send the same signal to the subprocess, the only difference is KeyBoardInterrupt is
        # manual, whereas timeout is a set variable
        for solver in self.solvers:
            with self.subTest("Test result of keyboard interrupt for solver", solver=solver.name):
                solver = solver(model=self.model)
                results = solver.run(timeout=1)
                self.assertFalse(results.to_array()[0][-1][0] == '5.0')
