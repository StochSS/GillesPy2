# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2024 GillesPy2 developers.

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
import math
import warnings
from unittest.mock import patch
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



class TestODEIntegrationFailure(unittest.TestCase):
    @staticmethod
    def create_transfer():
        model = gillespy2.Model(name='analytic_transfer')
        a = gillespy2.Species(name='A', initial_value=1, mode='continuous')
        b = gillespy2.Species(name='B', initial_value=0, mode='continuous')
        k = gillespy2.Parameter(name='k', expression=1)
        model.add_species([a, b])
        model.add_parameter(k)
        model.add_reaction(gillespy2.Reaction(name='transfer', reactants={a: 1},
                                            products={b: 1}, rate=k))
        model.timespan(np.linspace(0, 2, 3))
        return model

    def assert_transfer(self, result):
        trajectory = result[0]
        times = trajectory['time']
        np.testing.assert_allclose(trajectory['A'], np.exp(-times), rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(trajectory['B'], 1-np.exp(-times), rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(trajectory['A'] + trajectory['B'], 1, rtol=0, atol=1e-12)
        self.assertEqual(result[0].status, 'Success')

    def test_real_integrator_failure_modes_raise(self):
        for options, expected_code in [({'nsteps': 1}, -1),
                                       ({'rtol': 1e-30, 'atol': 1e-30}, -3)]:
            with self.subTest(options=options):
                solver = ODESolver(model=self.create_transfer())
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    with self.assertRaises(SimulationError) as caught:
                        solver.run(integrator='vode', integrator_options=options)
                self.assertEqual(solver.rc, -1)
                self.assertEqual(solver.result, [])
                self.assertIn('requested time 1.0', str(caught.exception.__cause__))
                self.assertIn('integrator return code ' + str(expected_code),
                              str(caught.exception.__cause__))

    def test_failure_after_valid_sample_does_not_publish_results(self):
        from gillespy2.solvers.numpy import ode_solver

        class FailedIntegrator:
            calls = []

            def __init__(self, function):
                self.t = 0

            def set_integrator(self, *args, **kwargs):
                return self

            def set_initial_value(self, y, t):
                self.t = t
                return self

            def set_f_params(self, *args):
                return self

            def integrate(self, requested):
                self.calls.append(requested)
                self.t = requested if len(self.calls) == 1 else 1.25
                return [math.exp(-self.t), 1-math.exp(-self.t)]

            def successful(self):
                return len(self.calls) == 1

            def get_return_code(self):
                return -5

        solver = ODESolver(model=self.create_transfer())
        with patch.object(ode_solver, 'ode', FailedIntegrator):
            with patch.object(ode_solver.Results, 'build_from_solver_results') as publish:
                with self.assertRaises(SimulationError) as caught:
                    solver.run()
                publish.assert_not_called()
        self.assertEqual(FailedIntegrator.calls, [1.0, 2.0])
        self.assertIn('requested time 2.0, reached 1.25', str(caught.exception.__cause__))
        self.assertEqual(solver.result, [])

    def test_failure_stops_at_first_output(self):
        from gillespy2.solvers.numpy import ode_solver
        backend_calls = []

        class FailedIntegrator:
            def __init__(self, function):
                self.t = 0.0

            def set_integrator(self, *args, **kwargs):
                return self

            def set_initial_value(self, y, t):
                return self

            def set_f_params(self, *args):
                return self

            def integrate(self, requested):
                backend_calls.append(requested)
                return [1, 0]

            def successful(self):
                return False

            def get_return_code(self):
                return -4

        solver = ODESolver(model=self.create_transfer())
        with patch.object(ode_solver, 'ode', FailedIntegrator):
            with self.assertRaises(SimulationError):
                solver.run()
        self.assertEqual(backend_calls, [1.0])
        self.assertEqual(solver.result, [])

    def test_successful_retry_after_failure(self):
        solver = ODESolver(model=self.create_transfer())
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with self.assertRaises(SimulationError):
                solver.run(integrator='vode', integrator_options={'nsteps': 1})
        self.assert_transfer(solver.run(integrator_options={'rtol': 1e-10, 'atol': 1e-12}))
        self.assertEqual(solver.rc, 0)
        self.assertFalse(hasattr(solver, 'has_raised_exception'))

    def test_failed_retry_does_not_return_prior_success(self):
        solver = ODESolver(model=self.create_transfer())
        self.assert_transfer(solver.run(integrator_options={'rtol': 1e-10, 'atol': 1e-12}))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            with self.assertRaises(SimulationError):
                solver.run(integrator='vode', integrator_options={'nsteps': 1})
        self.assertEqual(solver.result, [])

    def test_successful_analytic_solution_and_duplicate_trajectories(self):
        for integrator in ['vode', 'lsoda', 'dopri5', 'dop853']:
            with self.subTest(integrator=integrator):
                solver = ODESolver(model=self.create_transfer())
                result = solver.run(integrator=integrator, number_of_trajectories=2,
                                    integrator_options={'rtol': 1e-10, 'atol': 1e-12})
                self.assert_transfer(result)
                self.assertEqual(len(result), 2)
                np.testing.assert_array_equal(result[0]['A'], result[1]['A'])

    def test_resume_retains_valid_solution(self):
        model = self.create_transfer()
        solver = ODESolver(model=model)
        first = solver.run(integrator_options={'rtol': 1e-10, 'atol': 1e-12})
        resumed = solver.run(resume=first, t=4,
                             integrator_options={'rtol': 1e-10, 'atol': 1e-12})
        self.assert_transfer(resumed)
        self.assertEqual(resumed['time'][-1], 4)

if __name__ == '__main__':
    #unittest.main()
    TestBasicODESolver().test_stoch3()
