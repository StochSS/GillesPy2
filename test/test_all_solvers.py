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
from example_models import Example, Oregonator, MichaelisMenten
from gillespy2.core.results import Results, Trajectory
from gillespy2 import SSACSolver
from gillespy2 import ODESolver
from gillespy2 import NumPySSASolver
from gillespy2 import TauLeapingSolver
from gillespy2 import TauHybridSolver
from gillespy2 import ODECSolver
from gillespy2 import TauLeapingCSolver
from gillespy2 import TauHybridCSolver


class TestAllSolvers(unittest.TestCase):

    solvers = [
        SSACSolver,
        ODESolver,
        NumPySSASolver,
        TauLeapingSolver,
        TauHybridSolver,
        ODECSolver,
        TauLeapingCSolver,
        TauHybridCSolver,
    ]

    sbml_features = {
        "AssignmentRule": lambda model, variable:
            model.add_assignment_rule(gillespy2.AssignmentRule(variable=variable, formula="1/(t+1)")),
        "RateRule": lambda model, variable:
            model.add_rate_rule(gillespy2.RateRule(variable=variable, formula="2*t")),
        "Event": lambda model, variable:
            model.add_event(gillespy2.Event(
                trigger=gillespy2.EventTrigger(expression="t>1"),
                assignments=[gillespy2.EventAssignment(variable=variable, expression="100")]
            )),
        "FunctionDefinition": lambda model, variable:
            model.add_function_definition(
                gillespy2.FunctionDefinition(name="fn", function="variable", args=["variable"])),
    }

    # List of supported SBML features for each solver.
    # When a feature is implemented for a particular solver, add the feature to its list.
    solver_supported_sbml_features = {
        NumPySSASolver: [],
        TauLeapingSolver: [],
        ODESolver: [],
        TauHybridSolver: [
            "AssignmentRule",
            "RateRule",
            "Event",
            "FunctionDefinition",
        ],

        SSACSolver: [],
        ODECSolver: [],
        TauLeapingCSolver: [],
        TauHybridCSolver: [
            "RateRule",
            "Event",
        ],
    }

    model = Example()
    for sp in model.listOfSpecies.values():
        sp.mode = 'discrete'
    results = {}
    labeled_results = {}
    labeled_results_more_trajectories = {}

    for solver in solvers:
        labeled_results[solver] = model.run(solver=solver, number_of_trajectories=1,seed=1)
        labeled_results_more_trajectories[solver] = model.run(solver=solver, number_of_trajectories=2)

    def test_instantiated(self):
        for solver in self.solvers:
            self.model.run(solver=solver())

    def test_to_array(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver].to_array()[0], np.ndarray))

    def test_return_type_show_labels(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver], Results))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results[solver]['Sp'][0], np.float))

            self.assertTrue(isinstance(self.labeled_results[solver][0], Trajectory))

            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver], Results))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0], Trajectory))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver][0]['Sp'][0], np.float))

    def test_random_seed(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                same_results = self.model.run(solver=solver, seed=1)
                compare_results = self.model.run(solver=solver,seed=1)
                self.assertTrue(np.array_equal(same_results.to_array(), compare_results.to_array()))
                if solver.name in ["ODESolver", "ODECSolver"]: continue
                diff_results = self.model.run(solver=solver, seed=2)
                self.assertFalse(np.array_equal(diff_results.to_array(), same_results.to_array()))
    
    def test_random_seed_unnamed_reactions(self):
        model = self.model
        k2 = gillespy2.Parameter(name='k2', expression=1.0)
        model.add_parameter([k2])
        unnamed_rxn = gillespy2.Reaction(reactants={}, products={'Sp':1}, rate=k2)
        model.add_reaction(unnamed_rxn)
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                same_results = self.model.run(solver=solver, seed=1)
                compare_results = self.model.run(solver=solver,seed=1)
                self.assertTrue(np.array_equal(same_results.to_array(), compare_results.to_array()))
                if solver.name in ["ODESolver", "ODECSolver"]: continue
                diff_results = self.model.run(solver=solver, seed=2)
                self.assertFalse(np.array_equal(diff_results.to_array(), same_results.to_array()))

    def test_extraneous_args(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name), self.assertLogs(level='WARN'):
                model = Example()
                model.run(solver=solver, nonsense='ABC')

    def test_timeout(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name), self.assertLogs(level='WARN'):
                model = Oregonator()
                model.timespan(np.linspace(0, 1000000, 1001))
                model.run(solver=solver, timeout=0.1)

    def test_basic_solver_import(self):
        from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
        from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
        from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver
        model = MichaelisMenten()
        results1 = model.run(solver=BasicTauLeapingSolver)
        self.assertTrue(results1[0].solver_name == 'TauLeapingSolver')

        results2 = model.run(solver=BasicODESolver)
        self.assertTrue(results2[0].solver_name == 'ODESolver')

        results3 = model.run(solver=BasicTauHybridSolver)
        self.assertTrue(results3[0].solver_name == 'TauHybridSolver')

    @unittest.expectedFailure
    def test_sbml_feature_validation(self):
        class TestModel(gillespy2.Model):
            def __init__(self):
                gillespy2.Model.__init__(self, name="TestModel")
                self.add_species(gillespy2.Species(name="S", initial_value=0))
                self.timespan(np.linspace(0, 10, 11))

        all_features = set(self.sbml_features.keys())
        for solver in self.solvers:
            unsupported_features = all_features.difference(self.solver_supported_sbml_features.get(solver))
            with self.subTest(solver=solver.name):
                for sbml_feature_name in unsupported_features:
                    model = TestModel()
                    with self.subTest("Unsupported model features raise an error", sbml_feature=sbml_feature_name):
                        add_sbml_feature = self.sbml_features.get(sbml_feature_name)
                        add_sbml_feature(model, "S")
                        with self.assertRaises(gillespy2.ModelError):
                            solver.run(model=model)

                for sbml_feature_name in self.solver_supported_sbml_features.get(solver):
                    model = TestModel()
                    with self.subTest("Supported model features validate successfully", sbml_feature=sbml_feature_name):
                        add_sbml_feature = self.sbml_features.get(sbml_feature_name)
                        add_sbml_feature(model, "S")
                        solver.run(model=model)


if __name__ == '__main__':
    unittest.main()
