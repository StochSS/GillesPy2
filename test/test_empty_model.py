import unittest
import sys, os
sys.path.append(os.path.abspath(os.getcwd()))
import numpy, math, gillespy2
import matplotlib.pyplot as plt

class EmptyModel(gillespy2.Model):
    def __init__(self, parameter_values=None):
            gillespy2.Model.__init__(self)

class TestEmptyModel(unittest.TestCase):
    def setUp(self):
        self.model = EmptyModel()

    def test_model_creation(self):
        name = self.model.name
        self.assertEqual(name, "", msg="Unexpected value: {}".format(name))

if __name__ == '__main__':
    unittest.main()