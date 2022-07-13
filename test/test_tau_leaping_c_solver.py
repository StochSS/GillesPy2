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
import tempfile
from gillespy2.core.gillespyError import DirectoryError, SimulationError
from example_models import create_decay, create_decay_no_tspan
from gillespy2.solvers.cpp.tau_leaping_c_solver import TauLeapingCSolver

class TestTauLeapingCSolver(unittest.TestCase):
    def test_create(self):
        model = create_decay()
        solver = TauLeapingCSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = create_decay()
            solver = TauLeapingCSolver(model, temp.name)

    def test_run_example_precompiled(self):
        model = create_decay()
        solver = TauLeapingCSolver(model=model)
        results = model.run(solver=solver)

    def test_run_example(self):
        model = create_decay()
        solver = TauLeapingCSolver(model=model)
        results = model.run(solver=solver)

    def test_run_example__with_increment_only(self):
        model = create_decay_no_tspan()
        solver = TauLeapingCSolver(model=model)
        results = solver.run(increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = create_decay()
        solver = TauLeapingCSolver(model=model)
        results = solver.run()

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = create_decay()
            solver = TauLeapingCSolver(model=model)
            results = solver.run(increment=0.2)


if __name__ == '__main__':
    unittest.main()
