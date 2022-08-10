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
from example_models import create_decay
from gillespy2 import SSACSolver, ODECSolver, TauLeapingCSolver
from gillespy2 import TauHybridCSolver


class TestVariableSolvers(unittest.TestCase):
    model = create_decay()
    solverSSAC = SSACSolver(model, variable=True)
    solverODEC = ODECSolver(model, variable=True)
    solverTAUC = TauLeapingCSolver(model, variable=True)
    solverHYBC = TauHybridCSolver(model, variable=True)
    solverlist = [solverSSAC, solverODEC, solverTAUC, solverHYBC]

    def test_create(self):
        model = create_decay()
        solverSSAC = SSACSolver(model)
        solverODEC = ODECSolver(model)
        solverTAUC = TauLeapingCSolver(model)
        solverHYBC = TauHybridCSolver(model)

    def test_file_with_directory_name_exists(self):
        with self.assertRaises(DirectoryError):
            temp = tempfile.NamedTemporaryFile()
            model = create_decay()
            solverSSAC = SSACSolver(model, temp.name)
            solverODEC = ODECSolver(model, temp.name)
            solverTAUC = TauLeapingCSolver(model, temp.name)
            solverHYBC = TauHybridCSolver(model, temp.name)

    def test_run_example_precompiled(self):
        for solver in self.solverlist:
            results = self.model.run(solver=solver)

    def test_change_species(self):
        initial_value = self.model.listOfSpecies['Sp'].initial_value
        for solver in self.solverlist:
            results = self.model.run(solver=solver, variables={'Sp':3})
        with self.subTest(msg='Test changed species simulation'):
            self.assertEqual(results['Sp'][0], 3)
        with self.subTest(msg='Test changed species model integrity'):
            self.assertEqual(self.model.listOfSpecies['Sp'].initial_value, initial_value)

    def test_change_parameter(self):
        initial_expression = self.model.listOfParameters['k1'].expression
        for solver in self.solverlist:
            results = self.model.run(solver=solver, variables={'k1':0})
            with self.subTest(msg='Test changed parameter simulation'):
                self.assertEqual(results['Sp'][-1], results['Sp'][0])
            with self.subTest(msg='Test changed parameter model integrity'):
                self.assertEqual(self.model.listOfParameters['k1'].expression, initial_expression)

    def test_invalid_variable(self):
        with self.assertRaises(SimulationError):
            for solver in self.solverlist:
                results = self.model.run(solver=solver, variables={'foobar':0})


if __name__ == '__main__':
    unittest.main()
