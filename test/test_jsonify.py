import sys, unittest

sys.path.append("..")
from example_models import *

class TestJsonModels(unittest.TestCase):
    models = [
        Example,
        MichaelisMenten,
        VilarOscillator,
        Dimerization,
        Trichloroethylene,
        LacOperon,
        Schlogl,
        ToggleSwitch,
        Tyson2StateOscillator,
        Oregonator,
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
        from gillespy2.core.results import Results

        models = [
            Example,
            MichaelisMenten,
            Tyson2StateOscillator,
            Schlogl
        ]

        for model in models:
            target = model()
            results = target.run()

            r_json = results.to_json()
            r_from = Results.from_json(r_json)
            r_back = r_from.to_json()

            self.assertEqual(r_json, r_back)
