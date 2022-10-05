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
import numpy
from gillespy2.core import Model, Species, Reaction, Parameter, RateRule
from gillespy2.core.gillespyError import *


class TestSimpleModel(unittest.TestCase):
    def setUp(self):

        def create_simple_hybrid_model(parameter_values=None):
            model = Model(name="Simple_Hybrid_Model")

            # Species
            A = Species(name='A', initial_value=0)
            B = Species(name='B', initial_value=0)
            model.add_species([A, B])

            # Parameters
            k1 = Parameter(name='k1', expression=1)
            k2 = Parameter(name='k2', expression=10)
            model.add_parameter([k1, k2])

            # Rate Rule
            rate_rule = RateRule(name='Brate', variable='B', formula="cos(t)")
            model.add_rate_rule(rate_rule)

            # Reactions
            r1 = Reaction(name='r1', reactants={A: 1}, products={}, propensity_function="k1*B")
            r2 = Reaction(name='r2', reactants={}, products={B: 1}, rate=k2)
            model.add_reaction([r1, r2])

            model.timespan(numpy.linspace(0, 1, 11))
            return model

        self.model = create_simple_hybrid_model()

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
        self.assertIn('Name "{}" is unavailable. A species with that name exists.'.format(A.name), str(ex.exception))

    def test_addingMultipleSameSpecies_ThrowsError(self):
        A = Species(name='A', initial_value=0)
        B = Species(name='B', initial_value=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_species([A, B])
        self.assertIn('Name "{}" is unavailable. A species with that name exists.'.format(A.name), str(ex.exception))

    def test_addingSameParameter_ThrowsError(self):
        k1 = Parameter(name='k1', expression=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_parameter(k1)
        self.assertIn('Name "{}" is unavailable. A parameter with that name exists.'.format(k1.name), str(ex.exception))

    def test_addingMultipleSameParameter_ThrowsError(self):
        k1 = Parameter(name='k1', expression=0)
        k2 = Parameter(name='k2', expression=0)
        with self.assertRaises(ModelError) as ex:
            self.model.add_parameter([k1, k2])
        self.assertIn('Name "{}" is unavailable. A parameter with that name exists.'.format(k1.name), str(ex.exception))

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
        self.assertEqual(rate_rules['Brate'].variable.name, 'B', msg='Has incorrect species')
        self.assertEqual(rate_rules['Brate'].formula, 'cos(t)', msg='{0} has incorrect type'.format(rate_rules))

    def test_get_reaction(self):
        reaction = self.model.get_reaction('r1')
        self.assertIsInstance(reaction, Reaction, msg='{0} has incorrect type'.format(reaction))

    def test_model_has_reactions(self):
        reactions = self.model.get_all_reactions()
        self.assertIsInstance(reactions['r1'], Reaction, msg='{0} has incorrect type'.format(reactions))
        self.assertIsInstance(reactions['r2'], Reaction, msg='{0} has incorrect type'.format(reactions))

    def test_delete_all_reactions(self):
        self.model.delete_all_reactions()
        reactionList = self.model.get_all_reactions()
        self.assertNotIn('r1', reactionList)
        self.assertNotIn('r2', reactionList)

    def test_model_has_timespan_correct(self):
        timespan = self.model.tspan
        self.assertCountEqual(timespan, numpy.linspace(0, 1, 11), msg='Has incorrect timespan')


if __name__ == '__main__':
    unittest.main()
