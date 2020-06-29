import unittest
from example_models import Example
from gillespy2 import NumPySSASolver


class TestNumPySSASolver(unittest.TestCase):
    model = Example()
    

if __name__ == '__main__':
    unittest.main()
