import sys, unittest

sys.path.append("..")
from example_models import *

class TestJsonModels(unittest.TestCase):
    def test_jsonify_all_models(self):
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

        for model in models:
            self._test_model(model)


    def _test_model(self, model):
        print(f"Running: {model.__name__}")

        target = model()
        non_anon_json = target.to_json()
        non_anon_from_json = model.from_json(non_anon_json)
        non_anon_back_into_json = non_anon_from_json.to_json()

        self.assertEqual(non_anon_json, non_anon_back_into_json)
        print("\t[PASS] Non-anonymous json conversion.")

        target2 = model()
        table = target2.get_translation_table()

        anon_json = target2.to_json(table)
        anon_from_json = model.from_json(anon_json)
        anon_back_into_json = anon_from_json.to_json(table)

        self.assertEqual(anon_json, anon_back_into_json)
        print("\t[PASS] Anonymous json conversion.")
