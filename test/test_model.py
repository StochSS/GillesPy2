import unittest
from gillespy2.core import Model, Species, Reaction, Parameter
from gillespy2.core.gillespyError import *
import numpy as np


class TestModel(unittest.TestCase):

    def test_uniform_timespan(self):
        model = Model()
        model.timespan(np.linspace(0, 1, 100))
        with self.assertRaises(InvalidModelError):
            model.timespan(np.array([0, 0.1, 0.5]))

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

    def test_no_reaction_name(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=0)
        species2 = Species('B', initial_value=0)
        model.add_species([species1, species2])
        # add two reactions that has no name
        reaction1 = Reaction(reactants={species1: 1}, products={species2: 1}, rate=rate)
        reaction2 = Reaction(reactants={species2: 1}, products={species1: 1}, rate=rate)
        model.add_reaction([reaction1, reaction2])


    def test_int_type_mismatch(self):
        model = Model()
        y1 = np.int64(5)
        y2 = np.int32(5)
        species1 = Species('A', initial_value=y1)
        species2 = Species('B', initial_value=y2)
        model.add_species([species1, species2])

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

    def test_valid_initial_value_float(self):
        species = Species('A', initial_value=1.5, mode='continuous')

    def test_invalid_initial_value_float(self):
        with self.assertRaises(ValueError):
            species = Species('A', initial_value=1.5)

    def test_valid_initial_value_negative(self):
        species = Species('A', initial_value=-1, allow_negative_populations=True)

    def test_invalid_initial_value_negative(self):
        with self.assertRaises(ValueError):
            species = Species('A', initial_value=-1)

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
        results = model.run(number_of_trajectories=1, solver=NumPySSASolver, seed=1)
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
        with self.assertRaises(ParameterError):
            model.add_parameter(parameter)
        

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
        with self.assertRaises(Warning):
            model = Model(population = False, volume = 0.9)

    def test_model_init_custom_tspan(self):
        model = Model(tspan = np.linspace(0, 20, 401))

    def test_parameter_init_unspecified_expression(self):
        with self.assertRaises(TypeError):
            parameter = Parameter(name = 'parameter')

    def test_set_expression(self):
        parameter = Parameter(name = 'parameter', expression = 0.5)
        parameter.set_expression(1)

    def test_set_null_expression(self):
        parameter = Parameter(name = 'parameter', expression = 0.5)
        with self.assertRaises(TypeError):
            parameter.set_expression(None)

    def test_ode_propensity(self):
        model = Model()
        rate = Parameter(name='rate', expression=0.5)
        model.add_parameter(rate)
        species1 = Species('A', initial_value=10)
        species2 = Species('B', initial_value=10)
        model.add_species([species1, species2])
        r1 = Reaction(name='r1', reactants={'A':1}, products={}, rate=rate)
        r2 = Reaction(name='r2', reactants={'A':2}, products={'B':1}, rate=rate)
        r3 = Reaction(name='r3', reactants={'A':1, 'B':1}, products={}, rate=rate)
        r4 = Reaction(name='r4', reactants={'A':1}, products={}, propensity_function='t')
        model.add_reaction([r1, r2, r3, r4])
        self.assertEqual(model.listOfReactions['r1'].ode_propensity_function, 'rate*A')
        self.assertEqual(model.listOfReactions['r2'].ode_propensity_function, 'rate*A*A')
        self.assertEqual(model.listOfReactions['r3'].ode_propensity_function, 'rate*A*B')
        self.assertEqual(model.listOfReactions['r4'].ode_propensity_function, 't')

if __name__ == '__main__':
    unittest.main()
