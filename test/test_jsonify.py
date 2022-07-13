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

import sys
import unittest

sys.path.append("..")
from example_models import *
from gillespy2.core import Model, Reaction, Parameter, Species, Results
from gillespy2.core.jsonify import TranslationTable

class TestJsonModels(unittest.TestCase):
    models = [
        create_decay,
        create_michaelis_menten,
        create_vilar_oscillator,
        create_dimerization,
        create_trichloroethylene,
        create_toggle_switch,
        create_tyson_2_state_oscillator,
        create_oregonator,
    ]

    runnable_models = [
        create_decay,
        create_michaelis_menten,
        create_tyson_2_state_oscillator,
        create_schlogl
    ]

    def test_equality_of_named_models(self):
        """
        Test that a model can be converted to JSON and back.
        """

        for model in self.models:
            target = model()

            self.assertEqual(
                target, 
                Model.from_json(target.to_json())
            )

    def test_equality_of_anon_models(self):
        """
        Test that an anonymous model can be converted to JSON and back.
        """

        for model in self.models:
            target = model()
            anon_target = target.to_anon()

            self.assertEqual(
                anon_target, 
                Model.from_json(anon_target.to_json())
            )

    def test_equality_of_named_results(self):
        """
        Test that Model simulation results can be converted to JSON and back.
        """

        for model in self.runnable_models:
            target = model()
            results = target.run()

            results_from = Results.from_json(results.to_json())

            self.assertEquals(results.to_json(), results_from.to_json())

    def test_equality_of_named_results(self):
        """
        Test that anonymous Model simulation results can be converted to JSON and back.
        """

        for model in self.runnable_models:
            target = model()

            results = target.run()
            results._translation_table = target.get_translation_table()
            results = results.to_anon().to_json()

            self.assertEqual(results, Results.from_json(results).to_json())

    def test_model_hash_accuracy(self):
        """
        Test the accuracy of the JSON hash.
        """

        for model in self.models:
            # Create two two instances of the 'model' type.
            model_1 = model()
            model_2 = model()

            # Assert that the hash of the anonymized models are the same.
            self.assertEqual(
                model_1.to_anon().get_json_hash(), 
                model_2.to_anon().get_json_hash()
            )

            # Create a test class and change the variable insertion order.
            model_1 = model()
            model_1.var1 = "Hello"
            model_1.var2 = "world"
            model_1.var3 = [ "Hello world!" ]

            model_2 = model()
            model_2.var3 = [ "Hello world!" ]
            model_2.var2 = "world"
            model_2.var1 = "Hello"

            # Generate the first model's translation table.
            translation_table = model_1.get_translation_table()

            # Assert that the JSON of the anonymized models are still the same.
            self.assertEqual(
                model_1.to_anon().to_json(), 
                model_2.to_anon().to_json()
            )

            # Assert that the hash of the anonymized models are still the same.
            self.assertEqual(
                model_1.to_anon().get_json_hash(),
                model_2.to_anon().get_json_hash()
            )

            # Assert that the translation table is the same.
            self.assertEqual(
                model_1.get_translation_table().to_json(), 
                model_2.get_translation_table().to_json()
            )

            # Assert that model_1's JSON is equivalent to model_2 -> anon -> json -> object -> named -> json.
            self.assertEqual(
                model_1.to_json(), 
                Model.from_json(model_2.to_anon().to_json()).to_named().to_json()
            )

            # Assert that model_2's anon JSON hash is equivalent to model_2 -> anon -> json -> object -> json hash.
            self.assertEqual(
                model_1.to_anon().get_json_hash(), 
                Model.from_json(model_2.to_anon().to_json()).get_json_hash()
            )

    def test_named_to_anon_accuracy(self):
        for model in self.models:
            model_1 = model()

            # For each model, check to see if we can convert its table to and from json accurately.
            # For each model, generate its translation table and convert it to JSON.
            translation_table = model_1.get_translation_table()
            translation_table_json = translation_table.to_json()

            # Convert the JSON back into a TranslationTable object.
            translation_table_from_json = TranslationTable.from_json(translation_table_json)

            # Assert that the two tables are still identical.
            self.assertEqual(translation_table, translation_table_from_json)

            # Anonymize and convert model_1 to JSON.
            model_1 = model()
            model_1_json = model_1.to_anon().to_json()

            # Convert the JSON back into a Model object.
            model_2 = Model.from_json(model_1_json)

            # Assert that the anonymized model_1 and the new model_2 are identical.
            self.assertEquals(
                model_1.to_anon().to_json(),
                model_2.to_json()
            )

            # Convert the new model_2 to named.
            model_2 = model_2.to_named()

            # Assert that model_1 and model_2 are still the same.
            self.assertEquals(
                model_1.to_json(),
                model_2.to_json()
            )

    def test_model_hash_whitespace_accuracy(self):
        """ Test that differences in whitespace do not change the hash of a model. """
        model_no_whitespace = create_michaelis_menten()
        model_with_whitespace = create_michaelis_menten()

         
        X = Species(name="X", initial_value=int(0.65609071 * 300.0))
        Y = Species(name="Y", initial_value=int(0.85088331 * 300.0))

        model_no_whitespace.add_species([X, Y])
        model_with_whitespace.add_species([X, Y])

        # Up to this point the JSON hash of the two models should be the same.
        self.assertEquals(model_no_whitespace.get_json_hash(), model_with_whitespace.get_json_hash())

        # Add a custom reaction to both models, differing the amount of whitespace in the 
        # propensity functions.
        reaction_no_whitespace = Reaction(
            name="X production", 
            reactants={}, 
            products={X: 1},
            propensity_function="300*1.0/(1.0+(Y*Y/(300*300)))"
        )

        reaction_with_whitespace = Reaction(
            name="X production", 
            reactants={}, 
            products={X: 1},
            propensity_function="300      * 1.0 / (1.0 + (Y  *Y         /   (300  * 300)))"
        )

        model_no_whitespace.add_reaction(reaction_no_whitespace)
        model_with_whitespace.add_reaction(reaction_with_whitespace)
