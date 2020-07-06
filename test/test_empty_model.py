import unittest
from gillespy2.core.model import Model


class EmptyModel(Model):
    def __init__(self, parameter_values=None):
        Model.__init__(self)


class TestEmptyModel(unittest.TestCase):
    def setUp(self):
        self.model = EmptyModel()

    def test_model_creation(self):
        name = self.model.name
        self.assertEqual(name, "", msg="Unexpected value: {}".format(name))


if __name__ == '__main__':
    unittest.main()
