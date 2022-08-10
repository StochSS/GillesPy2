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
from example_models import create_robust_model, create_decay, create_decay_no_tspan
from gillespy2.core import Model, Species, Reaction, Parameter
from gillespy2.core.gillespyError import *
from gillespy2.core.model import export_StochSS
import tempfile
import os
import numpy as np


class TestModel(unittest.TestCase):

    def test_model_reflexivity(self):
        model = Model()
        assert model==model

    def test_model_inequality(self):
        model1 = Model()
        model2 = Model()
        param1 = Parameter('A', expression=0)
        param2 = Parameter('B', expression=1)
        model1.add_parameter(param1)
        model2.add_parameter(param2)
        assert model1 != model2

    def test_model_equality(self):
        model1 = Model()
        model2 = Model()
        param1 = Parameter('A', expression=0)
        model1.add_parameter(param1)
        model2.add_parameter(param1)
        assert model1 == model2

    def test_model_reordered_equality(self):
        model1 = Model()
        model2 = Model()
        param1 = Parameter('A', expression=0)
        param2 = Parameter('B', expression=1)
        model1.add_parameter(param1)
        model1.add_parameter(param2)
        model2.add_parameter(param2)
        model2.add_parameter(param1)
        assert model1 == model2

    def test_uniform_timespan(self):
        model = Model()
        model.timespan(np.linspace(0, 1, 100))
        with self.assertRaises(TimespanError):
            model.timespan(np.array([0, 0.1, 0.5]))
            model.run()

    def test_duplicate_parameter_names(self):
        model = Model()
        param1 = Parameter('A', expression=0)
        param2 = Parameter('A', expression=0)
        model.add_parameter(param1)
        with self.assertRaises(ModelError):
            model.add_parameter(param2)

    def test_duplicate_species_names(self):
        model = Model()
        species1 = Species('A', initial_value=0)
        species2 = Species('A', initial_value=0)
        model.add_species(species1)
        with self.assertRaises(ModelError):
            model.add_species(species2)

#    def test_int_type_mismatch(self):
#        model = Model()
#        y1 = np.int64(5)
#        y2 = np.int32(5)
#        species1 = Species('A', initial_value=y1)
#        species2 = Species('B', initial_value=y2)
#        model.add_species([species1, species2])

    def test_duplicate_reaction_name(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        reaction2 = Reaction(name="reaction1", reactants={species2: 1}, products={species1: 1}, rate=rate)
        model.add_reaction(reaction1)
        with self.assertRaises(ModelError):
            model.add_reaction(reaction2)

    def test_model_run_with_both_increment_and_timespan(self):
        model = create_decay()

        try:
            model.run(increment=4)

        except SimulationError as e:
            return

        self.fail(
            """
            Failed while testing Model.run() behavior when both `timespan` and `increment` are set.
            """
        )

    def test_reaction_invalid_reactant(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={'species1': 1}, products={species2: 1}, rate=rate)
        with self.assertRaises(ModelError):
            model.add_reaction(reaction1)

    def test_reaction_invalid_product(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={'species2': 1}, rate=rate)
        with self.assertRaises(ModelError):
            model.add_reaction(reaction1)
            
    def test_reaction_valid_reactant(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={'A': 1}, products={species2: 1}, rate=rate)
        model.add_reaction(reaction1)
        assert "reaction1" in model.listOfReactions

    def test_reaction_valid_product(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={'B': 1}, rate=rate)
        model.add_reaction(reaction1)
        assert "reaction1" in model.listOfReactions

    def test_add_reaction_dict(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        reactions = {
                name: Reaction(
                    name=name, 
                    reactants={species1: 1}, products={species2: 1}, rate=rate)
                for name in ["reaction1", "reaction2"]}
        with self.assertRaises(ModelError):
            model.add_reaction(reactions)
        
    def test_species_parameter_name_substrings(self):
        model = Model()
        rate = Parameter(name='rate', expression=1)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=100)
        species2 = Species('AA', initial_value=0)
        model.add_species([species1, species2])
        reaction1 = Reaction(name="reaction1", reactants={species1: 1}, products={species2: 1}, rate=rate)
        model.add_reaction(reaction1)
        number_points = 11
        model.timespan(np.linspace(0, 1, number_points))
        from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
        solver = NumPySSASolver(model=model)
        results = model.run(number_of_trajectories=1, solver=solver, seed=1)
        self.assertTrue(len(results['time']) == number_points)
        self.assertTrue(len(results[species1.name]) == number_points)
        self.assertTrue(len(results[species2.name]) == number_points)
        self.assertGreater(results[species1.name][0], results[species1.name][-1])
        self.assertLess(results[species2.name][0], results[species2.name][-1])
        self.assertEqual(np.sum(results[species1.name]) + np.sum(results[species2.name]), number_points * species1.initial_value)

    def test_problem_with_name(self):
        model = Model()
        numeric1 = Parameter(name = '123', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(numeric1)
        special1 = Parameter(name = 'parameter.', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special1)
        special2 = Parameter(name = 'parameter[', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special2)
        special3 = Parameter(name = 'parameter]', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special3)
        special4 = Parameter(name = 'parameter+', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special4)
        special5 = Parameter(name = 'parameter-', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special5)
        special6 = Parameter(name = 'parameter*', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special6)
        special7 = Parameter(name = 'parameter/', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special7)
        special8 = Parameter(name = 'parameter^', expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(special8)
        reserved1 = Parameter(name = "vol", expression = 0.5)
        with self.assertRaises(ModelError):
            model.add_parameter(reserved1)
        #Add this to name checks, python naming convention
        #numeric2 = Parameter(name = '1parameter', expression = 0.5)
        #with self.assertRaises(ModelError):
            #model.add_parameter(numeric2)

    def test_add_nonspecies_nonreaction_nonparameter(self):
        model = Model()
        species = 'nonspecies'
        with self.assertRaises(ModelError):
            model.add_species(species)
        reaction = 'nonreaction'
        with self.assertRaises(ModelError):
            model.add_reaction(reaction)
        parameter = 'nonparameter'
        with self.assertRaises(ModelError):
            model.add_parameter(parameter)
        
    def test_add_event(self):
        from gillespy2.core.events import Event, EventTrigger, EventAssignment
        model = Model()
        model.add_species(Species(name="A", initial_value=1, mode="discrete"))
        model.add_species(Species(name="B", initial_value=2, mode="discrete"))
        e1t = EventTrigger(expression="t>1")
        e1a1 = EventAssignment(variable="A", expression="3")
        e1a2 = EventAssignment(variable="B", expression="4")
        test_event=Event(name="e1", trigger=e1t, assignments=[e1a1, e1a2])
        model.add_event(test_event)
        passed_test = str(model)

    def test_run_nonsolver(self):
        model = Model()
        rate = Parameter(name = 'rate', expression = 0.5)
        model.add_parameter(rate)
        species1 = Species(name = 'A', initial_value = 0)
        species2 = Species(name = 'B', initial_value = 0)
        model.add_species(species1)
        model.add_species(species2)
        reaction = Reaction(name = 'reaction1', reactants={species1: 1}, products={species2: 1}, rate=rate)
        with self.assertRaises(SimulationError):
            results = model.run(number_of_trajectories = 1, solver = 'non_solver', seed = 1)

    def test_model_init_population_false_and_volume_warninag(self):
        with self.assertRaises(ModelError):
            model = Model(population = False, volume = 0.9)

    def test_model_init_custom_tspan(self):
        model = Model(tspan = np.linspace(0, 20, 401))

    def test_robust_model(self):
        try:
            model = create_robust_model()
            model.run()
        
        except ModelError as e:
            self.fail(f"Failed to instantiate the RobustModel: {e}")
        
        except SolverError as e:
            self.fail(f"Failed to run the RobustModel: {e}")

        except Exception as e:
            self.fail(f"An unknown exception occured while testing the RobustModel: {e}")

    def test_stochss_export(self):
        model = create_robust_model()
        tempdir = tempfile.mkdtemp()
        stochss_model_path = os.path.join(tempdir, "robust_model.mdl")
        try:
            export_StochSS(model, filename=stochss_model_path)
            self.assertTrue(os.path.exists(stochss_model_path))
        finally:
            os.unlink(stochss_model_path)
            os.rmdir(tempdir)

    def test_run_example__with_increment_only(self):
        model = create_decay_no_tspan()
        results = model.run(t=20, increment=0.2)

    def test_run_example__with_tspan_only(self):
        model = create_decay()
        results = model.run()

    def test_run_example__with_tspan_and_increment(self):
        with self.assertRaises(SimulationError):
            model = create_decay()
            results = model.run(t=20, increment=0.2)

if __name__ == '__main__':
    unittest.main()
