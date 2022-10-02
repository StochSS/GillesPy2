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
import re
import unittest

from gillespy2 import Model, Species, Parameter, Reaction
from gillespy2 import ReactionError

class TestReaction(unittest.TestCase):
    '''
    ################################################################################################
    Unit tests for gillespy2.Reaction.
    ################################################################################################
    '''
    def setUp(self):
        self.valid_ma_reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate="k1"
        )
        self.valid_cp_reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            propensity_function="k1 * A"
        )

    def test_constructor__mass_action(self):
        """ Test the Reaction constructor for a mass action reaction. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate="k1"
        )
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, {"A": 1})
        self.assertEqual(reaction.products, {"B": 1})
        self.assertEqual(reaction.marate, "k1")
        self.assertEqual(reaction.propensity_function, "(k1*A)")
        self.assertEqual(reaction.ode_propensity_function, "(k1*A)")

    def test_constructor__custom_propensity(self):
        """ Test the Reaction constructor for a custom propensity reaction. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            propensity_function="k1 * A * B"
        )
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, {"A": 1})
        self.assertEqual(reaction.products, {"B": 1})
        self.assertEqual(reaction.propensity_function, "((k1*A)*B)")
        self.assertEqual(reaction.ode_propensity_function, "((k1*A)*B)")

    def test_constructor__custom_ode_propensity(self):
        """ Test the Reaction constructor for a custom ode propensity reaction. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            ode_propensity_function="k1 * A * B"
        )
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, {"A": 1})
        self.assertEqual(reaction.products, {"B": 1})
        self.assertEqual(reaction.propensity_function, "((k1*A)*B)")
        self.assertEqual(reaction.ode_propensity_function, "((k1*A)*B)")

    def test_constructor__different_custom_propensity(self):
        """ Test the Reaction constructor for a custom propensity reaction with different propensities. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            propensity_function="k1 * A * B",
            ode_propensity_function="k1 * A * B / vol"
        )
        self.assertEqual(reaction.name, "test_reaction")
        self.assertEqual(reaction.reactants, {"A": 1})
        self.assertEqual(reaction.products, {"B": 1})
        self.assertEqual(reaction.propensity_function, "((k1*A)*B)")
        self.assertEqual(reaction.ode_propensity_function, "(((k1*A)*B)/vol)")

    def test_constructor__no_name(self):
        """ Test the Reaction constructor with no name provided. """
        reaction = Reaction(
            reactants={"A": 1}, products={"B": 1}, propensity_function="k1 * A * B"
        )
        self.assertIsNotNone(re.search("rxn.*", reaction.name))

    def test_constructor__name_is_none_or_empty(self):
        """ Test the Reaction constructor with None or empty string as name. """
        test_names = [None, ""]
        for test_name in test_names:
            with self.subTest(name=test_name):
                reaction = Reaction(
                    name=test_name, reactants={"A": 1}, products={"B": 1},
                    propensity_function="k1 * A * B"
                )
                self.assertIsNotNone(re.search("rxn.*", reaction.name))

    def test_constructor__invalid_name(self):
        """ Test the Reaction constructor with non-str name. """
        with self.assertRaises(ReactionError):
            Reaction(name=0, reactants={"A": 1}, products={"B": 1}, rate="k1")

    def test_constructor__no_reactants_or_products(self):
        """ Test the Reaction constructor with reactants and products not set. """
        with self.assertRaises(ReactionError):
            Reaction(name="test_reaction", rate="k1")

    def test_constructor__invalid_reactants_and_products(self):
        """ Test the Reaction constructor with reactants and products both set to None or empty. """
        test_reacs = [None, {}]
        test_prods = [None, {}]
        for test_reac in test_reacs:
            for test_prod in test_prods:
                with self.subTest(reactants=test_reac, products=test_prod):
                    with self.assertRaises(ReactionError):
                        Reaction(
                            name="test_reaction", reactants=test_reac,
                            products=test_prod, rate="k1"
                        )

    def test_constructor__reactants_is_none(self):
        """ Test the Reaction constructor with None as reactants. """
        reaction = Reaction(
            name="test_reaction", reactants=None, products={"B": 1}, rate="k2"
        )
        self.assertEqual(reaction.reactants, {})

    def test_constructor__reactants_keyed_by_species_obj(self):
        """ Test the Reaction constructor with reactants keyed by species objs. """
        test_species = Species(name="A", initial_value=0)
        reaction = Reaction(
            name="test_reaction", reactants={test_species: 1}, products={"B": 1},
            rate="k2"
        )
        self.assertEqual(reaction.reactants, {"A": 1})

    def test_constructor__reactants_invalid_keys(self):
        """ Test the Reaction constructor with reactants keyed with invalid keys. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={1: 1}, products={"B": 1},
                rate="k2"
            )

    def test_constructor__reactants_invalid_values(self):
        """ Test the Reaction constructor with reactants containing invalid stoichiometry. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": "1"}, products={"B": 1},
                rate="k2"
            )

    def test_constructor__invalid_reactants(self):
        """ Test the Reaction constructor with non-dict reactants. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants=[["A", 1]], products={"B": 1},
                rate="k2"
            )

    def test_constructor__products_is_none(self):
        """ Test the Reaction constructor with None as products. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products=None, rate="k2"
        )
        self.assertEqual(reaction.products, {})

    def test_constructor__products_keyed_by_species_obj(self):
        """ Test the Reaction constructor with products keyed by species objs. """
        test_species = Species(name="B", initial_value=0)
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={test_species: 1},
            rate="k2"
        )
        self.assertEqual(reaction.products, {"B": 1})

    def test_constructor__products_invalid_keys(self):
        """ Test the Reaction constructor with products keyed with invalid keys. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": 1}, products={1: 1},
                rate="k2"
            )

    def test_constructor__products_invalid_values(self):
        """ Test the Reaction constructor with products containing invalid stoichiometry. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": 1}, products={"B": "1"},
                rate="k2"
            )

    def test_constructor__invalid_products(self):
        """ Test the Reaction constructor with non-dict products. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": 1}, products=[["B", 1]],
                rate="k2"
            )

    def test_constructor__rate_and_propensity_functions_are_none(self):
        """ Test the Reaction constructor with rate and propensity functions set to None. """
        with self.assertRaises(ReactionError):
            Reaction(name="test_reaction", reactants={"A": 1}, products={"B": 1})

    def test_constructor__rate_and_propensity_function_are_not_none(self):
        """ Test the Reaction constructor with rate and propensity function set. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": 1}, products={"B": 1},
                rate="k1", propensity_function="A**2 + B**2"
            )

    def test_constructor__int_propensity_function(self):
        """ Test the Reaction constructor with an int propensity function. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            propensity_function=20
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "20")
        self.assertIsInstance(reaction.ode_propensity_function, str)
        self.assertEqual(reaction.ode_propensity_function, "20")

    def test_constructor__float_propensity_function(self):
        """ Test the Reaction constructor with a float propensity function. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            propensity_function=0.5
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "0.5")
        self.assertIsInstance(reaction.ode_propensity_function, str)
        self.assertEqual(reaction.ode_propensity_function, "0.5")

    def test_constructor__invalid_propensity_function(self):
        """ Test the Reaction constructor with a propensity function that is not of the proper type. """
        test_pfs = ["", ["k1 * A * B"]]
        for test_pf in test_pfs:
            with self.subTest(propensity_function=test_pf):
                with self.assertRaises(ReactionError):
                    Reaction(
                        name="test_reaction", reactants={"A": 1}, products={"B": 1},
                        propensity_function=test_pf
                    )

    def test_constructor__rate_and_ode_propensity_function_are_not_none(self):
        """ Test the Reaction constructor with rate and propensity function set. """
        with self.assertRaises(ReactionError):
            Reaction(
                name="test_reaction", reactants={"A": 1}, products={"B": 1},
                rate="k1", ode_propensity_function="A**2 + B**2"
            )

    def test_constructor__int_ode_propensity_function(self):
        """ Test the Reaction constructor with an int ode propensity function. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            ode_propensity_function=20
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "20")
        self.assertIsInstance(reaction.ode_propensity_function, str)
        self.assertEqual(reaction.ode_propensity_function, "20")

    def test_constructor__float_ode_propensity_function(self):
        """ Test the Reaction constructor with a float ode propensity function. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1},
            ode_propensity_function=0.5
        )
        self.assertIsInstance(reaction.propensity_function, str)
        self.assertEqual(reaction.propensity_function, "0.5")
        self.assertIsInstance(reaction.ode_propensity_function, str)
        self.assertEqual(reaction.ode_propensity_function, "0.5")

    def test_constructor__invalid_ode_propensity_function(self):
        """ Test the Reaction constructor with a ode propensity function that is not of the proper type. """
        test_opfs = ["", ["k1 * A * B"]]
        for test_opf in test_opfs:
            with self.subTest(ode_propensity_function=test_opf):
                with self.assertRaises(ReactionError):
                    Reaction(
                        name="test_reaction", reactants={"A": 1}, products={"B": 1},
                        ode_propensity_function=test_opf
                    )

    def test_constructor__parameter_rate(self):
        """ Test the Reaction constructor with an parameter object as rate. """
        k1 = Parameter(name="k1", expression="20")
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate=k1
        )
        self.assertIsInstance(reaction.marate, str)
        self.assertEqual(reaction.marate, "k1")

    def test_constructor__int_rate(self):
        """ Test the Reaction constructor with an int rate. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate=20
        )
        self.assertIsInstance(reaction.marate, str)
        self.assertEqual(reaction.marate, "20")

    def test_constructor__float_rate(self):
        """ Test the Reaction constructor with a float rate. """
        reaction = Reaction(
            name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate=0.5
        )
        self.assertIsInstance(reaction.marate, str)
        self.assertEqual(reaction.marate, "0.5")

    def test_constructor__rate_not_accepted_type(self):
        """ Test the Reaction constructor with a rate that is an invalid type. """
        test_rates = ["", ["k1"]]
        for test_rate in test_rates:
            with self.subTest(rate=test_rate):
                with self.assertRaises(ReactionError):
                    reaction = Reaction(
                        name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate=test_rate
                    )

    def test_constructor__annotation_invalid_type(self):
        """ Test the Reaction construct with an annotation of invalid type. """
        test_annotations = [5, 0.1, ["g"]]
        for test_annotation in test_annotations:
            with self.subTest(annotation=test_annotation):
                with self.assertRaises(ReactionError):
                    reaction = Reaction(
                        name="test_reaction", reactants={"A": 1}, products={"B": 1}, rate="k1",
                        annotation=test_annotation
                    )

    def test___str__(self):
        """ Test Reaction.__str__ method. """
        self.assertIsInstance(str(self.valid_ma_reaction), str)

    def test__create_mass_action__total_stoch_3(self):
        """ Test Reaction._create_mass_action total stochiometry > 2. """
        self.valid_ma_reaction.reactants = {"A": 1, "B": 2}
        self.valid_ma_reaction.products = {"C": 1}
        with self.assertRaises(ReactionError):
            self.valid_ma_reaction._create_mass_action()

    def test__create_mass_action__marate_type_as_string(self):
        """ Test Reaction._create_mass_action marate as string. """
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.products = {"C": 1}
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(k1*vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "k1")

    def test__create_mass_action__marate_type_as_int(self):
        """ Test Reaction._create_mass_action marate as int. """
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.products = {"C": 1}
        self.valid_ma_reaction.marate = 1
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(1*vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "1")

    def test__create_mass_action__marate_type_as_float(self):
        """ Test Reaction._create_mass_action marate as float. """
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.products = {"C": 1}
        self.valid_ma_reaction.marate = 0.5
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(0.5*vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "0.5")

    def test__create_mass_action__marate_type_as_parameter(self):
        """ Test Reaction._create_mass_action marate as parameter. """
        test_parameter = Parameter("k1", expression=0.1)
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.products = {"C": 1}
        self.valid_ma_reaction.marate = test_parameter
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(k1*vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "k1")

    def test__create_mass_action__X_to_Y(self):
        """ Test Reaction._create_mass_action X -> Y. """
        self.valid_ma_reaction.reactants = {"X": 1}
        self.valid_ma_reaction.products = {"Y": 1}
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(k1*X)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "(k1*X)")

    def test__create_mass_action__X_plus_Y_to_Z(self):
        """ Test Reaction._create_mass_action X + Y -> Z. """
        self.valid_ma_reaction.reactants = {"X": 1, "Y": 1}
        self.valid_ma_reaction.products = {"Z": 1}
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(((k1*X)*Y)/vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "((k1*X)*Y)")

    def test__create_mass_action__2X_to_Y(self):
        """ Test Reaction._create_mass_action 2X -> Y. """
        self.valid_ma_reaction.reactants = {"X": 2}
        self.valid_ma_reaction.products = {"Y": 1}
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(((k1*X)*(X-1))/vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "((k1*X)*X)")

    def test__create_mass_action__species_obj_reactant(self):
        """ Test Reaction._create_mass_action when reactants is keyed by species object. """
        test_species = Species(name="A", initial_value=1)
        self.valid_ma_reaction.reactants = {test_species: 1}
        self.valid_ma_reaction._create_mass_action()
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(k1*A)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "(k1*A)")

    def test__create_custom_propensity__exponent(self):
        """ Test Reaction._create_custom_propensity with a propensity containing exponents. """
        test_propensities = ["1/A^2", "1/A**2"]
        for test_propensity in test_propensities:
            with self.subTest(propensity_function=test_propensity):
                expr = self.valid_cp_reaction._create_custom_propensity(test_propensity)
                self.assertEqual(expr, "(1/pow(A,2))")

    def test__create_custom_propensity__ln(self):
        """ Test Reaction._create_custom_propensity with a propensity containing ln. """
        expr = self.valid_cp_reaction._create_custom_propensity("ln(A)")
        self.assertEqual(expr, "log(A)")

    def test__create_custom_propensity__e(self):
        """ Test Reaction._create_custom_propensity with a propensity containing e. """
        expr = self.valid_cp_reaction._create_custom_propensity("A*e")
        self.assertEqual(expr, "(A*2.718281828459045)")

    def test__create_custom_propensity__fake_e(self):
        """ Test Reaction._create_custom_propensity with a propensity containing e. """
        expr = self.valid_cp_reaction._create_custom_propensity("A*ex")
        self.assertEqual(expr, "(A*ex)")

    def test__create_custom_propensity__propensity_parsing(self):
        """ Test Reaction._create_custom_propensity for parsing accuracy. """
        test_propensities = [
            "5*x^2+e*b+6", "5*x**2+e*b+6", "1*alpha/2+5^beta", "1*alpha/2+5**beta",
            "2.78*x+3^(4*x)", "2.78*x+3**(4*x)", "-5*-x^2", "-5*-x**2",
            "(alpha/beta + delta**gamma)/(atlas-zeta)", "(alpha/beta + delta^gamma)/(atlas-zeta)"
        ]
        expected_results = [
            "(((5*pow(x,2))+(2.718281828459045*b))+6)", "(((5*pow(x,2))+(2.718281828459045*b))+6)",
            "(((1*alpha)/2)+pow(5,beta))", "(((1*alpha)/2)+pow(5,beta))", "((2.78*x)+pow(3,(4*x)))",
            "((2.78*x)+pow(3,(4*x)))", "((-5)*(-pow(x,2)))", "((-5)*(-pow(x,2)))",
            "(((alpha/beta)+pow(delta,gamma))/(atlas-zeta))", "(((alpha/beta)+pow(delta,gamma))/(atlas-zeta))"
        ]
        for i, test_propensity in enumerate(test_propensities):
            expected_result = expected_results[i]
            with self.subTest(propensity_function=test_propensity, expected_propensity_function=expected_result):
                expr = self.valid_cp_reaction._create_custom_propensity(test_propensity)
                self.assertEqual(expr, expected_result)

    def test_add_product__species_string(self):
        """ Test Reaction.add_product when species is string. """
        self.valid_ma_reaction.add_product("X", 1)
        self.assertIn("X", self.valid_ma_reaction.products)
        self.assertEqual(self.valid_ma_reaction.products["X"], 1)

    def test_add_product__species_object(self):
        """ Test Reaction.add_product when species is GillesPy2.Species. """
        test_species = Species(name="X", initial_value=1)
        self.valid_ma_reaction.add_product(test_species, 1)
        self.assertIn(test_species.name, self.valid_ma_reaction.products)
        self.assertEqual(self.valid_ma_reaction.products[test_species.name], 1)

    def test_add_product__invalid_species(self):
        """ Test Reaction.add_product with an invalid species. """
        test_species = [None, "", 5, 0.5, ["A"]]
        for test_spec in test_species:
            with self.subTest(species=test_spec):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.add_product(test_spec, 1)

    def test_add_product__invalid_stochiometry(self):
        """ Test Reaction.add_product with an invalid stochiometry. """
        test_stoichs = [None, "1", -5, 0, 0.5, [1]]
        for test_stoich in test_stoichs:
            with self.subTest(stoichiometry=test_stoich):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.add_product("X", test_stoich)

    def test_add_reactant__massaction_reaction_species_string(self):
        """ Test Reaction.add_reactant when reaction is mass-action and species is string. """
        self.valid_ma_reaction.add_reactant("X", 1)
        self.assertIn("X", self.valid_ma_reaction.reactants)
        self.assertEqual(self.valid_ma_reaction.reactants["X"], 1)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(((k1*A)*X)/vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "((k1*A)*X)")

    def test_add_reactant__customized_reaction_species_string(self):
        """ Test Reaction.add_reactant when reaction is customized and species is string. """
        self.valid_cp_reaction.add_reactant("X", 1)
        self.assertIn("X", self.valid_cp_reaction.reactants)
        self.assertEqual(self.valid_cp_reaction.reactants["X"], 1)
        self.assertEqual(self.valid_cp_reaction.propensity_function, "(k1*A)")
        self.assertEqual(self.valid_cp_reaction.ode_propensity_function, "(k1*A)")

    def test_add_reactant__massaction_reaction_species_object(self):
        """ Test Reaction.add_reactant when reaction is mass-action and species is GillesPy2.Species. """
        test_species = Species(name="X", initial_value=1)
        self.valid_ma_reaction.add_reactant(test_species, 1)
        self.assertIn(test_species.name, self.valid_ma_reaction.reactants)
        self.assertEqual(self.valid_ma_reaction.reactants[test_species.name], 1)
        self.assertEqual(self.valid_ma_reaction.propensity_function, f"(((k1*A)*{test_species.name})/vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, f"((k1*A)*{test_species.name})")

    def test_add_reactant__customized_reaction_species_object(self):
        """ Test Reaction.add_reactant when reaction is customized and species is GillesPy2.Species. """
        test_species = Species(name="X", initial_value=1)
        self.valid_cp_reaction.add_reactant(test_species, 1)
        self.assertIn(test_species.name, self.valid_cp_reaction.reactants)
        self.assertEqual(self.valid_cp_reaction.reactants[test_species.name], 1)
        self.assertEqual(self.valid_cp_reaction.propensity_function, "(k1*A)")
        self.assertEqual(self.valid_cp_reaction.ode_propensity_function, "(k1*A)")

    def test_add_reactant__massaction_reaction_stoich_0_to_1(self):
        """ Test Reaction.add_reactant when reaction is mass-action and stoichiometry goes from 0 to 1. """
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.add_reactant("X", 1)
        self.assertIn("X", self.valid_ma_reaction.reactants)
        self.assertEqual(self.valid_ma_reaction.reactants["X"], 1)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "(k1*X)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "(k1*X)")

    def test_add_reactant__massaction_reaction_stoich_0_to_2(self):
        """ Test Reaction.add_reactant when reaction is mass-action and stoichiometry goes from 0 to 2. """
        self.valid_ma_reaction.reactants = {}
        self.valid_ma_reaction.add_reactant("X", 2)
        self.assertIn("X", self.valid_ma_reaction.reactants)
        self.assertEqual(self.valid_ma_reaction.reactants["X"], 2)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "((((0.5*k1)*X)*(X-1))/vol)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "((k1*X)*X)")

    def test_add_reactant__invalid_species(self):
        """ Test Reaction.add_reactant with an invalid species. """
        test_species = [None, "", 5, 0.5, {}]
        for test_spec in test_species:
            with self.subTest(species=test_spec):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.add_reactant(test_spec, 1)

    def test_add_reactant__invalid_stochiometry(self):
        """ Test Reaction.add_reactant with an invalid stochiometry. """
        test_stoichs = [None, "1", -5, 0, 0.5, [1]]
        for test_stoich in test_stoichs:
            with self.subTest(stoichiometry=test_stoich):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.add_reactant("X", test_stoich)

    def test_from_json__massaction(self):
        """ Test Reaction.from_json with a mass-action reaction. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'A', 'value': 1}],
            'products': [{'key': 'B', 'value': 1}],
            'marate': 'k1',
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': True,
            'type': 'mass-action'
        }
        reaction = Reaction.from_json(test_json)
        self.assertIsInstance(reaction, Reaction)
        self.assertEqual(reaction.name, 'test_reaction')
        self.assertEqual(reaction.reactants, {'A': 1})
        self.assertEqual(reaction.products, {'B': 1})
        self.assertEqual(reaction.propensity_function, '(k1*A)')
        self.assertEqual(reaction.ode_propensity_function, '(k1*A)')
        self.assertEqual(reaction.marate, 'k1')
        self.assertIsNone(reaction.annotation)
        self.assertTrue(reaction.massaction)
        self.assertEqual(reaction.type, 'mass-action')

    def test_from_json__customized(self):
        """ Test Reaction.from_json with a customized reaction. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'A', 'value': 1}],
            'products': [{'key': 'B', 'value': 1}],
            'marate': None,
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': False,
            'type': 'customized'
        }
        reaction = Reaction.from_json(test_json)
        self.assertIsInstance(reaction, Reaction)
        self.assertEqual(reaction.name, 'test_reaction')
        self.assertEqual(reaction.reactants, {'A': 1})
        self.assertEqual(reaction.products, {'B': 1})
        self.assertEqual(reaction.propensity_function, '(k1*A)')
        self.assertEqual(reaction.ode_propensity_function, '(k1*A)')
        self.assertIsNone(reaction.marate)
        self.assertIsNone(reaction.annotation)
        self.assertFalse(reaction.massaction)
        self.assertEqual(reaction.type, 'customized')

    def test_from_json__massaction_json_changed(self):
        """ Test Reaction.from_json with a mass-action reaction and modified json. Ensure propensities are updated."""
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'B', 'value': 1}],
            'products': [{'key': 'A', 'value': 1}],
            'marate': 'k1',
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': True,
            'type': 'mass-action'
        }
        reaction = Reaction.from_json(test_json)
        self.assertEqual(reaction.name, 'test_reaction')
        self.assertEqual(reaction.reactants, {'B': 1})
        self.assertEqual(reaction.products, {'A': 1})
        self.assertEqual(reaction.propensity_function, '(k1*B)')
        self.assertEqual(reaction.ode_propensity_function, '(k1*B)')
        self.assertEqual(reaction.marate, 'k1')
        self.assertIsNone(reaction.annotation)
        self.assertTrue(reaction.massaction)
        self.assertEqual(reaction.type, 'mass-action')

    def test_from_json__customized_json_changed(self):
        """ Test Reaction.from_json with a customized reaction and modified json. Ensure propensities don't change. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'B', 'value': 1}],
            'products': [{'key': 'A', 'value': 1}],
            'marate': None,
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': False,
            'type': 'customized'
        }
        reaction = Reaction.from_json(test_json)
        self.assertEqual(reaction.name, 'test_reaction')
        self.assertEqual(reaction.reactants, {'B': 1})
        self.assertEqual(reaction.products, {'A': 1})
        self.assertEqual(reaction.propensity_function, '(k1*A)')
        self.assertEqual(reaction.ode_propensity_function, '(k1*A)')
        self.assertIsNone(reaction.marate)
        self.assertIsNone(reaction.annotation)
        self.assertFalse(reaction.massaction)
        self.assertEqual(reaction.type, 'customized')

    def test_from_json__invalid_json(self):
        """ Test Reaction.from_json with a reaction that fails validation. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'A', 'value': 1}],
            'products': [{'key': 'B', 'value': 1}],
            'marate': None,
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': True,
            'type': 'mass-action'
        }
        with self.assertRaises(ReactionError):
            reaction = Reaction.from_json(test_json)

    def test_set_annotation__invalid_annotation(self):
        """ Test Reaction.set_annotation with an invalid annotation. """
        test_annotations = [None, 5, 0.1, ["g"]]
        for test_annotation in test_annotations:
            with self.subTest(annotation=test_annotation):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.set_annotation(test_annotation)

    def test_set_propensities__propensity_function_only(self):
        """ Test Reaction.set_propensities with propensity_function only. """
        self.valid_ma_reaction.set_propensities(propensity_function="k1*A/2")
        expected_result = "((k1*A)/2)"
        self.assertEqual(self.valid_ma_reaction.propensity_function, expected_result)
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, expected_result)
        self.assertIsNone(self.valid_ma_reaction.marate)
        self.assertFalse(self.valid_ma_reaction.massaction)
        self.assertEqual(self.valid_ma_reaction.type, "customized")

    def test_set_propensities__ode_propensity_function_only(self):
        """ Test Reaction.set_propensities with ode_propensity_function only. """
        self.valid_ma_reaction.set_propensities(ode_propensity_function="k1*A/2")
        expected_result = "((k1*A)/2)"
        self.assertEqual(self.valid_ma_reaction.propensity_function, expected_result)
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, expected_result)
        self.assertIsNone(self.valid_ma_reaction.marate)
        self.assertFalse(self.valid_ma_reaction.massaction)
        self.assertEqual(self.valid_ma_reaction.type, "customized")

    def test_set_propensities__different_propensities(self):
        """ Test Reaction.set_propensities with different propensity functions. """
        self.valid_ma_reaction.set_propensities(
            propensity_function="k1*A/2", ode_propensity_function="k1*A/3"
        )
        self.assertEqual(self.valid_ma_reaction.propensity_function, "((k1*A)/2)")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "((k1*A)/3)")
        self.assertIsNone(self.valid_ma_reaction.marate)
        self.assertFalse(self.valid_ma_reaction.massaction)
        self.assertEqual(self.valid_ma_reaction.type, "customized")

    def test_set_propensities__int_propensity_function(self):
        """ Test Reaction.set_propensities with propensity_function of type int. """
        self.valid_ma_reaction.set_propensities(propensity_function=5)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "5")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "5")

    def test_set_propensities__float_propensity_function(self):
        """ Test Reaction.set_propensities with propensity_function of type float. """
        self.valid_ma_reaction.set_propensities(propensity_function=0.5)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "0.5")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "0.5")

    def test_set_propensities__int_ode_propensity_function(self):
        """ Test Reaction.set_propensities with ode_propensity_function of type int. """
        self.valid_ma_reaction.set_propensities(ode_propensity_function=5)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "5")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "5")

    def test_set_propensities__float_ode_propensity_function(self):
        """ Test Reaction.set_propensities with ode_propensity_function of type float. """
        self.valid_ma_reaction.set_propensities(propensity_function=0.5)
        self.assertEqual(self.valid_ma_reaction.propensity_function, "0.5")
        self.assertEqual(self.valid_ma_reaction.ode_propensity_function, "0.5")

    def test_set_propensities__invalid_propensity_function(self):
        """ Test Reaction.set_propensities with an invalid propensity_function. """
        test_pfs = ["", ["k1 * A * B"]]
        for test_pf in test_pfs:
            with self.subTest(propensity_function=test_pf):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.set_propensities(propensity_function=test_pf)

    def test_set_propensities__invalid_ode_propensity_function(self):
        """ Test Reaction.set_propensities with an invalid ode_propensity_function. """
        test_pfs = ["", ["k1 * A * B"]]
        for test_pf in test_pfs:
            with self.subTest(ode_propensity_function=test_pf):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.set_propensities(ode_propensity_function=test_pf)

    def test_set_rate(self):
        """ Test Reaction.set_rate. """
        self.valid_cp_reaction.set_rate("k2")
        self.assertEqual(self.valid_cp_reaction.propensity_function, "(k2*A)")
        self.assertEqual(self.valid_cp_reaction.ode_propensity_function, "(k2*A)")
        self.assertEqual(self.valid_cp_reaction.marate, "k2")
        self.assertTrue(self.valid_cp_reaction.massaction)
        self.assertEqual(self.valid_cp_reaction.type, "mass-action")

    def test_set_rate__parameter_rate(self):
        """ Test Reaction.set_rate with rate of type GillesPy2.Parameter. """
        test_parameter = Parameter(name="k2", expression="0.75")
        self.valid_cp_reaction.set_rate(test_parameter)
        self.assertEqual(self.valid_cp_reaction.marate, "k2")

    def test_set_rate__int_rate(self):
        """ Test Reaction.set_rate with rate of type int. """
        self.valid_cp_reaction.set_rate(5)
        self.assertEqual(self.valid_cp_reaction.marate, "5")

    def test_set_rate__float_rate(self):
        """ Test Reaction.set_rate with rate of type float. """
        self.valid_cp_reaction.set_rate(0.5)
        self.assertEqual(self.valid_cp_reaction.marate, "0.5")

    def test_set_rate__invalid_rate(self):
        """ Test Reaction.set_rate with an invalid rate. """
        test_rates = [None, "", ["k1"]]
        for test_rate in test_rates:
            with self.subTest(rate=test_rate):
                with self.assertRaises(ReactionError):
                    self.valid_cp_reaction.set_rate(test_rate)

    def test_to_dict__massaction(self):
        """ Test Reaction.to_dict with a mass-action reaction. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'A', 'value': 1}],
            'products': [{'key': 'B', 'value': 1}],
            'marate': 'k1',
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': True,
            'type': 'mass-action'
        }
        test_dict = self.valid_ma_reaction.to_dict()
        self.assertIsInstance(test_dict, dict)
        self.assertEqual(test_dict, test_json)

    def test_to_dict__customized(self):
        """ Test Reaction.to_dict with a customized reaction. """
        test_json = {
            'name': 'test_reaction',
            'reactants': [{'key': 'A', 'value': 1}],
            'products': [{'key': 'B', 'value': 1}],
            'marate': None,
            'annotation': None,
            'propensity_function': '(k1*A)',
            'ode_propensity_function': '(k1*A)',
            'massaction': False,
            'type': 'customized'
        }
        test_dict = self.valid_cp_reaction.to_dict()
        self.assertIsInstance(test_dict, dict)
        self.assertEqual(test_dict, test_json)

    def test_validate__mass_action(self):
        """ Test Reaction.validate for mass-action reactions. """
        test_types = ["mass-action", "customized"]
        test_massactions = [True, False]
        for test_type in test_types:
            for test_massaction in test_massactions:
                with self.subTest(reaction_type=test_type, massaction=test_massaction):
                    if test_type == "mass-action" and test_massaction:
                        continue
                    with self.assertRaises(ReactionError):
                        self.valid_ma_reaction.type = test_type
                        self.valid_ma_reaction.massaction = test_massaction
                        self.valid_ma_reaction.validate(coverage="all")

    def test_validate__customized(self):
        """ Test Reaction.validate for customized reactions. """
        test_types = ["mass-action", "customized"]
        test_massactions = [True, False]
        for test_type in test_types:
            for test_massaction in test_massactions:
                with self.subTest(reaction_type=test_type, massaction=test_massaction):
                    if test_type == "customized" and not test_massaction:
                        continue
                    with self.assertRaises(ReactionError):
                        self.valid_cp_reaction.type = test_type
                        self.valid_cp_reaction.massaction = test_massaction
                        self.valid_cp_reaction.validate(coverage="all")

    def test_validate__invalid_name(self):
        """ Test Reaction.validate with an invalid name. """
        test_names = [None, "", 0]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.name = test_name
                    self.valid_ma_reaction.validate(coverage="all")

    def test_valdiate__invalid_reactants_and_products(self):
        """ Test Reaction.validate with reactants and products both set to None or empty. """
        test_reacs = [None, {}]
        test_prods = [None, {}]
        for test_reac in test_reacs:
            for test_prod in test_prods:
                with self.subTest(reactants=test_reac, products=test_prod):
                    with self.assertRaises(ReactionError):
                        self.valid_ma_reaction.reactants = test_reac
                        self.valid_ma_reaction.products = test_prod
                        self.valid_ma_reaction.validate(coverage="all")

    def test_validate__reactants_invalid_keys(self):
        """ Test Reaction.validate with reactants keyed with invalid keys. """
        with self.assertRaises(ReactionError):
            self.valid_ma_reaction.reactants = {1: 1}
            self.valid_ma_reaction.validate(coverage="all")

    def test_validate__reactants_invalid_values(self):
        """ Test Reaction.validate with reactants containing invalid stoichiometry. """
        with self.assertRaises(ReactionError):
            self.valid_ma_reaction.reactants = {"A": "1"}
            self.valid_ma_reaction.validate(coverage="all")

    def test_validate__invalid_reactants(self):
        """ Test Reaction.validate with non-dict reactants. """
        test_reactants = [None, ["A", 1]]
        for test_reactant in test_reactants:
            with self.subTest(reactants=test_reactant):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.reactants = test_reactant
                    self.valid_ma_reaction.validate(coverage="all")

    def test_validate__prodects_invalid_keys(self):
        """ Test Reaction.validate with products keyed with invalid keys. """
        with self.assertRaises(ReactionError):
            self.valid_ma_reaction.products = {1: 1}
            self.valid_ma_reaction.validate(coverage="all")

    def test_validate__products_invalid_values(self):
        """ Test Reaction.validate with products containing invalid stoichiometry. """
        with self.assertRaises(ReactionError):
            self.valid_ma_reaction.products = {"B": "1"}
            self.valid_ma_reaction.validate(coverage="all")

    def test_validate__invalid_products(self):
        """ Test Reaction.validate with non-dict products. """
        test_products = [None, ["B", 1]]
        for test_product in test_products:
            with self.subTest(products=test_product):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.products = test_product
                    self.valid_ma_reaction.validate(coverage="all")

    def test_validate__invalid_propensity_function(self):
        """ Test Reaction.validate with an invalid propensity function. """
        test_pfs = [None, 20, 0.5, "", ["k1 * A * B"]]
        for test_pf in test_pfs:
            with self.subTest(propensity_function=test_pf):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.propensity_function = test_pf
                    self.valid_ma_reaction.validate(coverage="all")

    def test_validate__invalid_ode_propensity_function(self):
        """ Test Reaction.validate with an invalid ode propensity function. """
        test_opfs = [None, 20, 0.5, "", ["k1 * A * B"]]
        for test_opf in test_opfs:
            with self.subTest(ode_propensity_function=test_opf):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.ode_propensity_function = test_opf
                    self.valid_ma_reaction.validate(coverage="all")

    def test_validate__invalid_rate(self):
        """ Test Reaction.validate with a rate that is an invalid type. """
        test_rates = [20, 0.5, "", ["k1"]]
        for test_rate in test_rates:
            with self.subTest(rate=test_rate):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.marate = test_rate
                    self.valid_ma_reaction.validate(coverage="all")

    def test_validate__annotation_invalid_type(self):
        """ Test Reaction.validate with an annotation of invalid type. """
        test_annotations = [5, 0.1, ["g"]]
        for test_annotation in test_annotations:
            with self.subTest(annotation=test_annotation):
                with self.assertRaises(ReactionError):
                    self.valid_ma_reaction.annotation = test_annotation
                    self.valid_ma_reaction.validate(coverage="all")

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        import time
        from datetime import datetime
        start = time.time()
        self.valid_ma_reaction.validate(coverage="all")
        tic = datetime.utcfromtimestamp(time.time() - start)
        print(f"Total time to run validate on a mass-action reaction: {tic.strftime('%M mins %S secs %f msecs')}")
        start = time.time()
        self.valid_cp_reaction.validate(coverage="all")
        tic = datetime.utcfromtimestamp(time.time() - start)
        print(f"Total time to run validate on a customized reaction: {tic.strftime('%M mins %S secs %f msecs')}")
