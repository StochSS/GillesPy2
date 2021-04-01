import sys, unittest

sys.path.append("..")
from example_models import *

class TestJsonModels(unittest.TestCase):
    def __init__(self):
        import inspect

        models = list((x for x in inspect.importlib.import_module("example_models") if x.isclass))
        for model in models:
            self.test_model(model)