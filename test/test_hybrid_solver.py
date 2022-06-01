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
import gillespy2
from gillespy2.core.gillespyError import *
from example_models import create_decay, create_decay_no_tspan, create_multi_firing_event
from gillespy2 import TauHybridSolver


class TestBasicTauHybridSolver(unittest.TestCase):

    def test_null_timeout(self):
        model = create_decay()
        species = gillespy2.Species('test_species', initial_value=1)
        model.add_species([species])
        solver = TauHybridSolver(model=model)
        results = solver.run()
        self.assertTrue(len(results) > 0)

    def test_add_rate_rule(self):
        model = create_decay()
        species = gillespy2.Species('test_species', initial_value=1)
        rule = gillespy2.RateRule(name='rr1',formula='test_species+1',variable='test_species')
        model.add_species([species])
        model.add_rate_rule([rule])
        results = model.run()
        valid_solvers = ('TauHybridSolver', 'TauHybridCSolver')
        self.assertIn(results[0].solver_name, valid_solvers)

    def test_add_rate_rule_dict(self):
        model = create_decay()
        species2 = gillespy2.Species('test_species2', initial_value=2, mode='continuous')
        species3 = gillespy2.Species('test_species3', initial_value=3, mode='continuous')
        rule2 = gillespy2.RateRule('rule2', species2, 'cos(t)')
        rule3 = gillespy2.RateRule(name='rule3', variable=species3, formula='sin(t)')
        rate_rule_dict = {'rule2': rule2, 'rule3': rule3}
        model.add_species([species2, species3])
        with self.assertRaises(ModelError):
            model.add_rate_rule(rate_rule_dict)

    def test_add_bad_species_rate_rule_dict(self):
        model = create_decay()
        rule = gillespy2.RateRule(formula='sin(t)')
        with self.assertRaises(ModelError):
            model.add_rate_rule(rule)

    def test_add_assignment_rule(self):
        model = create_decay()
        species = gillespy2.Species('test_species4', initial_value=1)
        rule = gillespy2.AssignmentRule(name='ar1', variable=species.name, formula='2')
        model.add_species([species])
        model.add_assignment_rule([rule])
        results = model.run()
        self.assertEquals(results[species.name][0], 2) 
        self.assertEquals(results[species.name][-1], 2)
        self.assertEqual(results[0].solver_name,'TauHybridSolver')

    def test_add_function_definition(self):
        model = create_decay()
        funcdef = gillespy2.FunctionDefinition(name='fun', function='Sp+1')
        model.add_function_definition(funcdef)
        results = model.run()
        self.assertEqual(results[0].solver_name,'TauHybridSolver')

    def test_add_continuous_species_dependent_event(self):
        model = create_decay()
        model.listOfSpecies['Sp'].mode = 'continuous'
        eventTrig = gillespy2.EventTrigger(expression='Sp <= 90', initial_value=True, )
        event1 = gillespy2.Event(name='event1', trigger=eventTrig)
        ea1 = gillespy2.EventAssignment(variable='Sp', expression='1000')
        ea2 = gillespy2.EventAssignment(variable='k1', expression='0')
        event1.add_assignment([ea1, ea2])
        model.add_event(event1)
        results = model.run()
        valid_solvers = ('TauHybridSolver', 'TauHybridCSolver')
        self.assertIn(results[0].solver_name, valid_solvers)
        self.assertEqual(results['Sp'][-1], 1000)

    def test_add_stochastic_species_dependent_event(self):
        model = create_decay()
        model.listOfSpecies['Sp'].mode = 'discrete'
        eventTrig = gillespy2.EventTrigger(expression='Sp <= 90', initial_value=True, )
        event1 = gillespy2.Event(name='event1', trigger=eventTrig)
        ea1 = gillespy2.EventAssignment(variable='Sp', expression='1000')
        ea2 = gillespy2.EventAssignment(variable='k1', expression='0')
        event1.add_assignment([ea1, ea2])
        model.add_event(event1)
        results = model.run()
        self.assertEqual(results['Sp'][-1], 1000)
        
    def test_add_continuous_time_dependent_event(self):
        model = create_decay()
        model.listOfSpecies['Sp'].mode = 'continuous'
        eventTrig = gillespy2.EventTrigger(expression='t >= 10', initial_value=True, )
        event1 = gillespy2.Event(name='event1', trigger=eventTrig)
        ea1 = gillespy2.EventAssignment(variable='Sp', expression='1000')
        ea2 = gillespy2.EventAssignment(variable='k1', expression='0')
        event1.add_assignment([ea1, ea2])
        model.add_event(event1)
        results = model.run()
        self.assertEqual(results['Sp'][-1], 1000)
        
    def test_add_stochastic_time_dependent_event(self):
        model = create_decay()
        model.listOfSpecies['Sp'].mode = 'discrete'
        eventTrig = gillespy2.EventTrigger(expression='t >= 10', initial_value=True, )
        event1 = gillespy2.Event(name='event1', trigger=eventTrig)
        ea1 = gillespy2.EventAssignment(variable='Sp', expression='1000')
        ea2 = gillespy2.EventAssignment(variable='k1', expression='0')
        event1.add_assignment([ea1, ea2])
        model.add_event(event1)
        results = model.run()
        self.assertEqual(results['Sp'][-1], 1000)
        
    def test_add_param_event(self):
        def create_event_test_model(parameter_values=None):
            model = gillespy2.Model(name='Event Test Model')
            model.add_species([gillespy2.Species(name='S', initial_value=0)])
            model.add_parameter(gillespy2.Parameter(name='event_tracker',
                                                    expression=99))
            model.add_parameter(gillespy2.Parameter(name='event_tracker2',
                                                    expression=0))
            model.add_reaction(gillespy2.Reaction(name='r1', products={'S':1},
                                                rate=model.listOfParameters['event_tracker2']))
            eventTrig1 = gillespy2.EventTrigger(expression='t>=2')
            event1 = gillespy2.Event(name='event1', trigger=eventTrig1)
            event1.add_assignment(gillespy2.EventAssignment(
                                    variable='event_tracker', expression='t'))
            eventTrig2 = gillespy2.EventTrigger(expression='t >= event_tracker + 2')
            event2 = gillespy2.Event(name='event2', trigger=eventTrig2)
            event2.add_assignment(gillespy2.EventAssignment(
                                    variable='event_tracker2', expression='t'))
            model.add_event([event1, event2])
            return model

        model = create_event_test_model()
        results = model.run(increment=0.05, t=20)
        self.assertGreater(results['S'][-1], 0)

    def test_math_name_overlap(self):
        model = create_decay()
        gamma = gillespy2.Species('gamma',initial_value=2, mode='continuous')
        model.add_species([gamma])
        k2 = gillespy2.Parameter(name='k2', expression=1)
        model.add_parameter([k2])
        gamma_react = gillespy2.Reaction(name='gamma_react', reactants={'gamma': 1}, products={}, rate=k2)
        model.add_reaction([gamma_react])
        solver = TauHybridSolver(model=model)
        model.run(solver=solver)

    def test_add_bad_expression_rate_rule_dict(self):
        model = create_decay()
        species2 = gillespy2.Species('test_species2', initial_value=2, mode='continuous')
        rule = gillespy2.RateRule(variable=species2, formula='')
        with self.assertRaises(ModelError):
            model.add_rate_rule(rule)

    def test_ensure_hybrid_dynamic_species(self):
        model = create_decay()
        species1 = gillespy2.Species('test_species1', initial_value=1,mode='dynamic')
        model.add_species(species1)
        results = model.run()
        valid_solvers = ('TauHybridSolver', 'TauHybridCSolver')
        self.assertIn(results[0].solver_name, valid_solvers)

    def test_ensure_hybrid_continuous_species(self):
        model = create_decay()
        species1 = gillespy2.Species('test_species1', initial_value=1,mode='continuous')
        model.add_species(species1)
        results = model.run()
        valid_solvers = ('TauHybridSolver', 'TauHybridCSolver')
        self.assertIn(results[0].solver_name, valid_solvers)

    def test_ensure_continuous_dynamic_timeout_warning(self):
        model = create_decay()
        species1 = gillespy2.Species('test_species1', initial_value=1, mode='dynamic')
        model.add_species(species1)
        with self.assertLogs('GillesPy2', level='WARN'):
            solver = TauHybridSolver(model=model)
            results = model.run(solver=solver, timeout=1)

    def test_run_example__with_increment_only(self):
        model = create_decay_no_tspan()
        solver = TauHybridSolver(model=model)
        results = solver.run(increment=0.2, t=20)

    def test_run_example__with_tspan_only(self):
        model = create_decay()
        solver = TauHybridSolver(model=model)
        results = solver.run()

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = create_decay()
            solver = TauHybridSolver(model=model)
            results = solver.run(increment=0.2)


class TestAllHybridSolvers(unittest.TestCase):
    from gillespy2.solvers import TauHybridCSolver

    solvers = [
        TauHybridSolver,
        TauHybridCSolver,
    ]

    def test_continuous_state_values(self):
        """
        Continuous values should be evaluate appropriately, without being truncated/casted to an integer.
        """
        
        def create_truncated_state_model(parameter_values=None):
            model = gillespy2.Model(name="TruncatedStateModel")
            S1 = gillespy2.Species(name="S1", initial_value=0, mode="discrete")
            rate = gillespy2.Species(name="rate", initial_value=0.9999, mode="continuous")
            model.add_species([S1, rate])
            model.add_rate_rule(gillespy2.RateRule(variable="rate", formula="-1/((t+0.9999)**2)"))
            model.add_reaction(
                # Because S1 is a "discrete" species, our reaction will be marked "stochastic."
                gillespy2.Reaction(products={S1: 1}, propensity_function="10.0*rate")
            )
            return model

        model = create_truncated_state_model()
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                solver = solver(model=model)
                result = model.run(solver=solver, seed=1, increment=0.05, t=20)
                self.assertGreater(result["S1"][-1], 0.0,
                                   "Reaction never fired; indicates that continuous species is being truncated")

    def test_multi_firing_event(self):
        model = create_multi_firing_event()
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                solver = solver(model=model)
                res = model.run(solver=solver, seed=1)
                self.assertNotEqual(res['Sp'][45], 0)
                self.assertEqual(res['Sp'][75], 0)
                self.assertNotEqual(res['Sp'][96], 0)
                self.assertEqual(res['Sp'][120], 0)
                self.assertNotEqual(res['Sp'][144], 0)
                self.assertEqual(res['Sp'][165], 0)

if __name__ == '__main__':
    unittest.main()
