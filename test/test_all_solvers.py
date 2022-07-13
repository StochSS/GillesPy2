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
from example_models import create_decay, create_oregonator, create_michaelis_menten
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

    @classmethod
    def setUpClass(cls):
        cls.solvers = [
            SSACSolver,
            ODESolver,
            NumPySSASolver,
            TauLeapingSolver,
            TauHybridSolver,
            ODECSolver,
            TauLeapingCSolver,
            TauHybridCSolver,
        ]

        cls.sbml_features = {
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
        cls.solver_supported_sbml_features = {
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

        cls.model = create_decay()

        cls.results = {}
        cls.labeled_results = {}
        cls.labeled_results_more_trajectories = {}
    
    def test_extraneous_args(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name), self.assertLogs('GillesPy2', level='WARN'):
                model = create_decay()
                solver = solver(model=model)
                model.run(solver=solver, nonsense='ABC')

    def test_instantiated(self):
        for sp in self.model.listOfSpecies.values():
            sp.mode = 'discrete'

        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                solver = solver(model=self.model)
                if "ODE" in solver.name:
                    self.labeled_results[solver.name] = self.model.run(solver=solver, number_of_trajectories=1)
                else:
                    self.labeled_results[solver.name] = self.model.run(solver=solver, number_of_trajectories=1, seed=1)
                self.labeled_results_more_trajectories[solver.name] = self.model.run(solver=solver, number_of_trajectories=2)

    def test_random_seed(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                solver = solver(model=self.model)
                if "ODE" in solver.name:
                    same_results = self.model.run(solver=solver)
                    compare_results = self.model.run(solver=solver)
                else:
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
                solver = solver(model=self.model)
                if "ODE" in solver.name:
                    same_results = self.model.run(solver=solver)
                    compare_results = self.model.run(solver=solver)
                else:
                    same_results = self.model.run(solver=solver, seed=1)
                    compare_results = self.model.run(solver=solver,seed=1)
                self.assertTrue(np.array_equal(same_results.to_array(), compare_results.to_array()))
                if solver.name in ["ODESolver", "ODECSolver"]: continue
                diff_results = self.model.run(solver=solver, seed=2)
                self.assertFalse(np.array_equal(diff_results.to_array(), same_results.to_array()))

    def test_return_type_show_labels(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver.name], Results))
            self.assertTrue(isinstance(self.labeled_results[solver.name]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results[solver.name]['Sp'][0], np.float))

            self.assertTrue(isinstance(self.labeled_results[solver.name][0], Trajectory))

            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver.name], Results))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver.name][0], Trajectory))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver.name][0]['Sp'], np.ndarray))
            self.assertTrue(isinstance(self.labeled_results_more_trajectories[solver.name][0]['Sp'][0], np.float))

    def test_sbml_feature_validation(self):
        def create_test_model(parameter_values=None):
            model = gillespy2.Model(name="TestModel")
            model.add_species(gillespy2.Species(name="S", initial_value=0))
            model.timespan(np.linspace(0, 10, 11))
            return model

        all_features = set(self.sbml_features.keys())
        for solver in self.solvers:
            unsupported_features = all_features.difference(self.solver_supported_sbml_features.get(solver))
            with self.subTest(solver=solver.name):
                for sbml_feature_name in unsupported_features:
                    model = create_test_model()
                    with self.subTest("Unsupported model features raise an error", sbml_feature=sbml_feature_name):
                        add_sbml_feature = self.sbml_features.get(sbml_feature_name)
                        add_sbml_feature(model, "S")
                        with self.assertRaises(gillespy2.ModelError):
                            solver.validate_sbml_features(model=model)

                for sbml_feature_name in self.solver_supported_sbml_features.get(solver):
                    model = create_test_model()
                    with self.subTest("Supported model features validate successfully", sbml_feature=sbml_feature_name):
                        add_sbml_feature = self.sbml_features.get(sbml_feature_name)
                        add_sbml_feature(model, "S")
                        solver.validate_sbml_features(model=model)

    def test_timeout(self):
        for solver in self.solvers:
            with self.subTest(solver=solver.name), self.assertLogs('GillesPy2', level='WARN'):
                model = create_oregonator()
                model.timespan(np.linspace(0, 1000000, 1001))
                solver = solver(model=model)
                model.run(solver=solver, timeout=0.1)

    def test_to_array(self):
        for solver in self.solvers:
            self.assertTrue(isinstance(self.labeled_results[solver.name].to_array()[0], np.ndarray))


if __name__ == '__main__':
    unittest.main()
