import unittest
import numpy
from gillespy2.core import Model, Species, Reaction, Parameter, RateRule
from gillespy2.core.gillespyError import *


class SimpleHybridModel(Model):
    def __init__(self, parameter_values=None):
        Model.__init__(self, name="Simple_Hybrid_Model")

        # Species
        A = Species(name='A', initial_value=0)
        B = Species(name='B', initial_value=0)
        self.add_species([A, B])

        # Parameters
        k1 = Parameter(name='k1', expression=1)
        k2 = Parameter(name='k2', expression=10)
        self.add_parameter([k1, k2])

        # Rate Rule
        rate_rule = RateRule(name='Brate', variable='B', formula="cos(t)")
        self.add_rate_rule(rate_rule)

        # Reactions
        r1 = Reaction(name='r1', reactants={A: 1}, products={}, propensity_function="k1*B")
        r2 = Reaction(name='r2', reactants={}, products={B: 1}, rate=k2)
        self.add_reaction([r1, r2])

        self.timespan(numpy.linspace(0, 1, 11))


class TestSimpleModel(unittest.TestCase):
    def setUp(self):
        self.model = SimpleHybridModel()

#    def test_this_should_fail(self):
#        name = self.model.name
#        self.assertEqual(name, "NOT Simple_Hybrid_Model", msg="Unexpected value: {}".format(name))

    def test_model_creation(self):
        name = self.model.name
        self.assertEqual(name, "Simple_Hybrid_Model", msg="Unexpected value: {}".format(name))

    def test_addingSameSpecies_ThrowsError(self):
        A = Species(name='A', initial_value=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_species(A)
        self.assertEqual(str(ex.exception), 'Name "{}" is unavailable. A species with that name exists.'.format(A.name))

    def test_addingMultipleSameSpecies_ThrowsError(self):
        A = Species(name='A', initial_value=0)
        B = Species(name='B', initial_value=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_species([A, B])
        self.assertEqual(str(ex.exception), 'Name "{}" is unavailable. A species with that name exists.'.format(A.name))

    def test_addingSameParameter_ThrowsError(self):
        k1 = Parameter(name='k1', expression=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_parameter(k1)
        self.assertEqual(str(ex.exception), 'Name "{}" is unavailable. A parameter with that name exists.'.format(k1.name))

    def test_addingMultipleSameParameter_ThrowsError(self):
        k1 = Parameter(name='k1', expression=0)
        k2 = Parameter(name='k2', expression=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_parameter([k1, k2])
        self.assertEqual(str(ex.exception), 'Name "{}" is unavailable. A parameter with that name exists.'.format(k1.name))

    def test_delete_species(self):
        self.model.delete_species('A')
        speciesList = self.model.get_all_species()
        self.assertNotIn('A', speciesList)

    def test_delete_all_species(self):
        self.model.delete_all_species()
        speciesList = self.model.get_all_species()
        self.assertNotIn('A', speciesList)
        self.assertNotIn('B', speciesList)

    def test_set_units_concentration(self):
        self.model.set_units('concentration')
        units = self.model.units
        self.assertEqual(units, 'concentration')

    def test_set_units_population(self):
        self.model.set_units('population')
        units = self.model.units
        self.assertEqual(units, 'population')

    def test_setFakeUnits_ThrowsError(self):
        with self.assertRaises(ModelError) as ex:
            self.model.set_units('nonsense')
        self.assertEqual(str(ex.exception), "units must be either concentration or population (case insensitive)")

    def test_serialization(self):
        doc = self.model.serialize()
        self.assertIsInstance(doc, str)
        self.assertIn('<Reaction>', doc)

    def test_model_has_species(self):
        species = self.model.get_species('A')
        self.assertIsInstance(species, Species, msg='{0} has incorrect type'.format(species))

    def test_add_single_species(self):
        C = Species(name='C', initial_value=0)
        self.model.add_species(C)
        species = self.model.get_species('C')
        self.assertIsInstance(species, Species, msg='{0} has incorrect type'.format(species))

    def test_model_has_species_list(self):
        speciesList = self.model.get_all_species()
        species = speciesList['A']
        self.assertIsInstance(species, Species, msg='{0} has incorrect type'.format(species))

    def test_get_parameter(self):
        parameter = self.model.get_parameter('k1')
        self.assertEqual('1', parameter.expression)
        self.assertIsInstance(parameter, Parameter, msg='{0} has incorrect type'.format(parameter))

    def test_getFakeParameter_ThrowsError(self):
        p_name = 'fake'
        with self.assertRaises(ModelError) as ex:
            parameter = self.model.get_parameter(p_name)
        self.assertEqual(str(ex.exception), "No parameter named " + p_name)

    def test_set_parameter(self):
        self.model.set_parameter('k1', '100')
        parameter = self.model.get_parameter('k1')
        self.assertEqual('100', parameter.expression)
        self.assertIsInstance(parameter, Parameter, msg='{0} has incorrect type'.format(parameter))

    def test_delete_parameter(self):
        self.model.delete_parameter('k1')
        parameterList = self.model.get_all_parameters()
        self.assertNotIn('k1', parameterList)

    def test_delete_all_parameters(self):
        self.model.delete_all_parameters()
        parameterList = self.model.get_all_parameters()
        self.assertNotIn('k1', parameterList)
        self.assertNotIn('k2', parameterList)

    def test_model_has_parameters(self):
        parameters = self.model.get_all_parameters()
        self.assertIsInstance(parameters['k1'], Parameter, msg='{0} has incorrect type'.format(parameters))
        self.assertIsInstance(parameters['k2'], Parameter, msg='{0} has incorrect type'.format(parameters))

    def test_model_parameters_correct(self):
        parameters = self.model.get_all_parameters()
        self.assertEqual(parameters['k1'].expression, '1', msg='Has incorrect expression')
        self.assertEqual(parameters['k2'].expression, '10', msg='Has incorrect expression')

    def test_model_has_rate_rules(self):
        rate_rules = self.model.listOfRateRules
        print(rate_rules)
        self.assertEqual(rate_rules['B'].variable, 'B', msg='Has incorrect species')
        self.assertEqual(rate_rules['B'].formula, 'cos(t)', msg='{0} has incorrect type'.format(rate_rules))

    def test_get_reaction(self):
        reaction = self.model.get_reaction('r1')
        self.assertIsInstance(reaction, Reaction, msg='{0} has incorrect type'.format(reaction))

    def test_model_has_reactions(self):
        reactions = self.model.get_all_reactions()
        self.assertIsInstance(reactions['r1'], Reaction, msg='{0} has incorrect type'.format(reactions))
        self.assertIsInstance(reactions['r2'], Reaction, msg='{0} has incorrect type'.format(reactions))

    def test_delete_reaction(self):
        self.model.delete_reaction('r1')
        reactionList = self.model.get_all_reactions()
        self.assertNotIn('r1', reactionList)

    def test_delete_all_reactions(self):
        self.model.delete_all_reactions()
        reactionList = self.model.get_all_reactions()
        self.assertNotIn('r1', reactionList)
        self.assertNotIn('r2', reactionList)

    def test_model_has_reactions_correct(self):
        reactions = self.model.get_all_reactions()

        species_A = self.model.get_species('A')
        species_B = self.model.get_species('B')

        reactants_r1 = reactions['r1'].reactants
        products_r2 = reactions['r2'].products

        species_r1 = list(reactants_r1)[0]
        species_r2 = list(products_r2)[0]

        # Check r1 name & propensity function is set
        self.assertEqual(reactions['r1'].name, 'r1', msg='Has incorrect expression')
        self.assertEqual(reactions['r1'].propensity_function, 'k1*B', msg='Has incorrect expression')

        # Check r1 reactants are set
        self.assertEqual(reactants_r1[species_A], 1, msg='Has incorrect number of reactants')
        self.assertIsInstance(species_r1, Species, msg='Has incorrect type')

        # Check r2 products are set
        self.assertEqual(products_r2[species_B], 1, msg='Has incorrect number of products')
        self.assertIsInstance(species_r2, Species, msg='Has incorrect type')

        # Check r2 name & rate is set
        self.assertEqual(reactions['r2'].name, 'r2', msg='Has incorrect expression')
        self.assertEqual(reactions['r2'].marate.expression, '10', msg='Has incorrect expression')

    def test_model_has_timespan_correct(self):
        timespan = self.model.tspan
        self.assertCountEqual(timespan, numpy.linspace(0, 1, 11), msg='Has incorrect timespan')


if __name__ == '__main__':
    unittest.main()
