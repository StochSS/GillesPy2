import unittest

def capital_case(x):
    return x.capitalize()

class TestSimpleModel(unittest.TestCase):
    def test_capital_case(self):
        self.assertEqual('Semaphore', capital_case('semaphore'))

if __name__ == '__main__':
    unittest.main()