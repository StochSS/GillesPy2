from gillespy2.core.reaction import Reaction
from gillespy2.core import parameter
from gillespy2.core.parameter import Parameter
import sys, unittest

sys.path.append("..")
from test.example_models import *
from gillespy2.core.results import Results

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
        for model in self.models:
            target = model()
            non_anon_json = target.to_json()
            non_anon_from_json = model.from_json(non_anon_json)
            non_anon_back_into_json = non_anon_from_json.to_json()

            self.assertEqual(non_anon_json, non_anon_back_into_json)

    def test_anon_model_norun(self):
        for model in self.models:
            target = model()
            table = target.get_translation_table()

            anon_json = target.to_json(table)
            anon_from_json = model.from_json(anon_json)
            anon_back_into_json = anon_from_json.to_json(table)

            self.assertEqual(anon_json, anon_back_into_json)

    def test_non_anon_model_run(self):
        for model in self.runnable_models:
            target = model()
            results = target.run()

            r_json = results.to_json()
            r_from = Results.from_json(r_json)
            r_back = r_from.to_json()

            self.assertEqual(r_json, r_back)

    def test_anon_model_runs(self):
        for model in self.runnable_models:
            target = model()
            translation_table = target.get_translation_table()
            results = target.run()

            r_json = results.to_json(translation_table)
            r_from = Results.from_json(r_json)
            r_back = r_from.to_json(translation_table)

            self.assertEqual(r_json, r_back)

    def test_model_hash(self):
        for model in self.models:
            # Simple test to see if two identical models will return the same hash.
            model_1 = model()
            model_2 = model()

            self.assertEqual(model_1.get_json_hash(model_1.get_translation_table()), model_2.get_json_hash(model_2.get_translation_table()))

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

            # A bit overkill, but it's good to check to ensure that both the JSON and hash output match.
            self.assertEqual(model_1.to_json(translation_table), model_2.to_json(translation_table))
            self.assertEqual(model_1.get_json_hash(translation_table), model_2.get_json_hash(translation_table))

    def test_model_hash_chaos(self):
        for model in self.models:
            model_1 = model()
            model_2 = model()

            import random
            from gillespy2.core import Parameter, Species, Reaction

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

            translation_table = model_1.get_translation_table()
            self.assertEqual(model_1.get_json_hash(translation_table), model_2.get_json_hash(translation_table))

            # The translation_table of model_2 should still make model_1 and model_2 equivalent.
            translation_table = model_2.get_translation_table()
            self.assertEqual(model_1.get_json_hash(translation_table), model_2.get_json_hash(translation_table))