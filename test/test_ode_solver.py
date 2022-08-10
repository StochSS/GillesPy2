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
import sys
sys.path.append("..")
import numpy as np
import gillespy2
from example_models import create_decay, create_decay_no_tspan
from gillespy2 import ODESolver
from gillespy2.core.gillespyError import SimulationError


class TestBasicODESolver(unittest.TestCase):
    number_of_trajectories = 10

    def test_run_example(self):
        model = create_decay()
        solver = ODESolver(model=model)
        for i in [1, self.number_of_trajectories]:
            for label in [True, False]:
                with self.subTest(number_of_trajectories=i, show_labels=label):
                    if i > 1:
                        results = model.run(solver=solver, number_of_trajectories=i)
                        self.assertEqual(len(results), i)
                    else:
                        results = model.run(solver=solver, number_of_trajectories=i)

                    if not label:
                        results = results.to_array()
                    if i > 1 or not label:
                        for result in results[1:]:
                            if label:
                                self.assertEqual(results[0], result)
                            else:
                                self.assertTrue(np.array_equal(results[0], result))

    def test_stoich2(self):
        def create_stoch_test_1(parameter_values=None):
            model = gillespy2.Model(name='StochTest1')
            A = gillespy2.Species(name='A', initial_value=10)
            B = gillespy2.Species(name='B', initial_value=0)
            model.add_species([A, B])
            k = gillespy2.Parameter(name='k', expression=10)
            model.add_parameter([k])
            r = gillespy2.Reaction(name='r', reactants={A: 2}, products={B:1}, rate=k)
            model.add_reaction([r])
            model.timespan(np.linspace(0, 100, 101))
            return model
    
        model = create_stoch_test_1()
        solver = ODESolver(model=model)
        result = model.run(solver=solver)
        self.assertAlmostEqual(result['B'][-1], 5, places=3)

    def test_run_example__with_increment_only(self):
        model = create_decay_no_tspan()
        solver = ODESolver(model=model)
        results = solver.run(increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = create_decay()
        solver = ODESolver(model=model)
        results = solver.run()

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = create_decay()
            solver = ODESolver(model=model)
            results = solver.run(increment=0.2)


    def test_stoch3(self):
        def create_stoch_test_1(parameter_values=None):
            model = gillespy2.Model(name='StochTest1')
            A = gillespy2.Species(name='A', initial_value=10)
            B = gillespy2.Species(name='B', initial_value=0)
            model.add_species([A, B])
            k = gillespy2.Parameter(name='k', expression=10)
            model.add_parameter([k])
            r = gillespy2.Reaction(name='r', reactants={A: 1}, products={B:1},
                propensity_function="k*A/vol") # testing if 'vol' is a pre-set variable
            model.add_reaction([r])
            model.timespan(np.linspace(0, 100, 101))
            return model

        model = create_stoch_test_1()
        solver = ODESolver(model=model)
        result = model.run(solver=solver)
        sys.stderr.write(f"\ntest_shoch3(): B={result['B'][-1]}\n\n")
        self.assertGreater(result['B'][-1], 5)

if __name__ == '__main__':
    #unittest.main()
    TestBasicODESolver().test_stoch3()
