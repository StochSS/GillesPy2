"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import numpy as np
import gillespy2
from example_models import Example, ExampleNoTspan
from gillespy2 import ODESolver
from gillespy2.core.gillespyError import SimulationError


class TestBasicODESolver(unittest.TestCase):
    number_of_trajectories = 10

    def test_run_example(self):
        model = Example()
        for i in [1, self.number_of_trajectories]:
            for label in [True, False]:
                with self.subTest(number_of_trajectories=i, show_labels=label):
                    if i > 1:
                        with self.assertLogs(level='WARN'):
                            results = model.run(solver=ODESolver, show_labels=label, number_of_trajectories=i)
                        self.assertEqual(len(results), i)
                    else:
                        results = model.run(solver=ODESolver, show_labels=label, number_of_trajectories=i)

                    if i > 1 or not label:
                        for result in results[1:]:
                            if label:
                                self.assertEqual(results[0], result)

                            else:
                                self.assertTrue(np.array_equal(results[0], result))

    def test_stoich2(self):
        class StoichTestModel(gillespy2.Model):
            def __init__(self, parameter_values=None):
                gillespy2.Model.__init__(self, name='StochTest1')
                A = gillespy2.Species(name='A', initial_value=10)
                B = gillespy2.Species(name='B', initial_value=0)
                self.add_species([A, B])
                k = gillespy2.Parameter(name='k', expression=10)
                self.add_parameter([k])
                r = gillespy2.Reaction(name='r', reactants={A: 2}, products={B:1}, rate=k)
                self.add_reaction([r])
                self.timespan(np.linspace(0, 100, 101))
        model = StoichTestModel()
        result = model.run(solver=ODESolver)
        self.assertAlmostEqual(result['B'][-1], 5, places=3)

    def test_run_example__with_increment_only(self):
        model = ExampleNoTspan()
        results = ODESolver.run(model=model, increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = Example()
        results = ODESolver.run(model=model)

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = Example()
            results = ODESolver.run(model=model, increment=0.2)



if __name__ == '__main__':
    unittest.main()
