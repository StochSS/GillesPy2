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
import subprocess
from shutil import which
from example_models import create_michaelis_menten, create_oregonator
from gillespy2.core.results import Results, Trajectory
from gillespy2.core import Species
from gillespy2 import SSACSolver
from gillespy2 import ODESolver
from gillespy2 import NumPySSASolver
from gillespy2 import TauLeapingSolver
from gillespy2.core import gillespyError
import subprocess
import signal
import time
import os


class TestPauseResume(unittest.TestCase):
    solvers = [SSACSolver, ODESolver, NumPySSASolver, TauLeapingSolver]

    model = create_michaelis_menten()
    for sp in model.listOfSpecies.values():
        sp.mode = 'discrete'
    labeled_results = {}
    labeled_results_more_trajectories = {}


    for solver in solvers:
        solver = solver(model=model)
        labeled_results[solver.name] = model.run(solver=solver)

    def test_altered_model_failure(self):
        model = create_michaelis_menten()
        for solver_class in self.solvers:
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
            self.labeled_results[solver.name] = model.run(solver=solver,
                                                     resume=self.labeled_results[solver.name], t=150)
        for solver in self.solvers:
            self.assertEqual(int(self.labeled_results[solver.name][0]['time'][-1]),150)

    def test_time_fail(self):
        model = self.model
        for solver in self.solvers:
            with self.assertRaises((gillespyError.ExecutionError, gillespyError.SimulationError)):
                solver = solver(model=model)
                self.labeled_results = model.run(solver=solver, resume=self.labeled_results[solver.name], t=1)

    def test_pause(self):
        py_path = which('python3')
        if py_path is None:
            py_path = which('python')
        model_path = os.path.join(os.path.dirname(__file__), 'pause_model.py')
        args = [[py_path, model_path, 'NumPySSASolver'],
                [py_path, model_path, 'TauLeapingSolver'],
                [py_path, model_path, 'ODESolver']]
        for arg in args:
            if os.name == 'nt':
                p = subprocess.Popen(arg, stdout=subprocess.PIPE, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)
                time.sleep(2)
                p.send_signal(signal.CTRL_BREAK_EVENT)
            else:
                p = subprocess.Popen(arg, start_new_session=True, stdout=subprocess.PIPE)
                time.sleep(2)
                os.kill(p.pid, signal.SIGINT)
            out, err = p.communicate()
            # End time for Oregonator is 5. If indexing into a numpy array using the form:
            # results[0][-1][0] (where .run(show_labels=False), this index being the last time in the index
            # One would get an output of "5.0" before converting it to an int. Hence, assert time != 5.0 rather than 5.
            self.assertFalse(out.decode('utf-8').rstrip() == '5.0')

        solvers = [SSACSolver]
        # For the C solvers, timeouts behave identical to a keyboard interrupt, and would return the same data, if
        # one was to KeyBoardInterrupt, at the same time a timeout ended. This is because timeouts in C solvers
        # and KeyBoardInterrupt send the same signal to the subprocess, the only difference is KeyBoardInterrupt is
        # manual, whereas timeout is a set variable
        for solver in solvers:
            model = create_oregonator()
            solver = solver(model=model)
            results = model.run(solver=solver, timeout=1)
            self.assertFalse(results.to_array()[0][-1][0] == '5.0')



