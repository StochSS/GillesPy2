import sys
import unittest

sys.path.append("..")
from example_models import *
from gillespy2.core import Model, Reaction, Parameter, Species, Results
from gillespy2.core.jsonify import TranslationTable

class TestJsonModels(unittest.TestCase):
    models = [
        Example,
        MichaelisMenten,
        VilarOscillator,
        Dimerization,
        Trichloroethylene,
        ToggleSwitch,
        Tyson2StateOscillator,
        Oregonator,
    ]

    runnable_models = [
        Example,
        MichaelisMenten,
        Tyson2StateOscillator,
        Schlogl
    ]

    def test_non_anon_model_norun(self):
        """
        Test that a model can be converted to JSON and back.
        """

        for model in self.models:
            target = model()

            self.assertEqual(
                target, 
                model.from_json(target.to_json())
            )

    def test_anon_model_norun(self):
        """
        Test that an anonymous model can be converted to JSON and back.
        """

        for model in self.models:
            target = model()
            anon_target = target.to_anon()

            self.assertEqual(
                anon_target, 
                model.from_json(anon_target.to_json())
            )

    def test_non_anon_model_run(self):
        """
        Test that Model simulation results can be converted to JSON and back.
        """

        for model in self.runnable_models:
            target = model()
            results = target.run()

            results_from = Results.from_json(results.to_json())

            self.assertEquals(results.to_json(), results_from.to_json())

    def test_anon_model_runs(self):
        """
        Test that anonymous Model simulation results can be converted to JSON and back.
        """

        for model in self.runnable_models:
            target = model()

            results = target.run()
            results._translation_table = target.get_translation_table()
            results = results.to_anon().to_json()

            self.assertEqual(results, Results.from_json(results).to_json())

    def test_model_hash(self):
        """
        Test the accuracy of the JSON hash.
        """

        for model in self.models:
            # Simple test to see if two identical models will return the same hash.
            model_1 = model()
            model_2 = model()

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

            translation_table = model_1.get_translation_table()

            # Test to ensure that the JSON output and JSON hash are equal.
            self.assertEqual(
                model_1.to_anon().to_json(), 
                model_2.to_anon().to_json()
            )

            self.assertEqual(
                model_1.to_anon().get_json_hash(),
                model_2.to_anon().get_json_hash()
            )

            # Test to ensure that the translation table for both models is identical.
            self.assertEqual(
                model_1.get_translation_table().to_json(), 
                model_2.get_translation_table().to_json()
            )

            # Test to ensure that model_1's JSON is equivalent to model_2 -> anon -> json -> object -> named -> json.
            self.assertEqual(
                model_1.to_json(), 
                Model.from_json(model_2.to_anon().to_json()).to_named().to_json()
            )

            # Test to ensure that model_2's anon JSON hash is equivalent to model_2 -> anon -> json -> object -> json hash.
            self.assertEqual(
                model_1.to_anon().get_json_hash(), 
                Model.from_json(model_2.to_anon().to_json()).get_json_hash()
            )

    def test_model_hash_chaos(self):
        import random

        for model in self.models:
            model_1 = model()
            model_2 = model()

            # Generate lists of random parameters, species, and reactions. Each will then be added to model_1 and model_2 in random order.
            parameters = [Parameter(name=bytes(random.sample(range(97, 123), 10)).decode(), expression=random.randint(0, 10)) for x in range(5)]
            model_1.add_parameter(random.sample(parameters, len(parameters)))
            model_2.add_parameter(random.sample(parameters, len(parameters)))

            species = [Species(name=bytes(random.sample(range(97, 123), 10)).decode(), initial_value=random.randint(0, 10)) for x in range(5)]
            model_1.add_species(random.sample(species, len(species)))
            model_2.add_species(random.sample(species, len(species)))

            reactions = [Reaction(
                name=bytes(random.sample(range(97, 123), 10)).decode(), 
                reactants={ species[random.randint(0, len(species) - 1)].name: random.randint(0, 5) },
                products={ species[random.randint(0, len(species) - 1)].name: random.randint(0, 5) },
                propensity_function=parameters[random.randint(0, len(parameters) - 1)].name) for x in range(5)]

            model_1.add_reaction(random.sample(reactions, len(reactions)))
            model_2.add_reaction(random.sample(reactions, len(reactions)))

            # At this point, model_1 and model_2 contain the same data, but it was entered in a different order.
            # The json hash function should ensure that they are still equivalent.
            self.assertEqual(model_1.to_anon().get_json_hash(), model_2.to_anon().get_json_hash())

            # The translation table for model_1 and model_2 should also be the same.
            self.assertEqual(model_1.get_translation_table().to_json(), model_2.get_translation_table().to_json())

    def test_anon_conversion(self):
        for model in self.models:
            model_1 = model()

            # For each model, check to see if we can convert its table to and from json accurately.
            translation_table = model_1.get_translation_table()
            translation_table_from_json = TranslationTable.from_json(translation_table.to_json())

            self.assertEqual(translation_table.to_json(), translation_table_from_json.to_json())

            # For each model, ensure that the anonymized version is still equivalent to the original when converted back.
            model_1 = model()
            model_2 = model.from_json(model_1.to_anon().to_json()).to_named()

            self.assertEqual(model_1.to_json(), model_2.to_json())