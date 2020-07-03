import unittest
import sys

class TestInit(unittest.TestCase):
    
    def test_sys_init(self):
        x = sys.version_info
        sys.version_info = (2,7)

        with self.assertRaises(Exception):
            Example().run()

        sys.version_info = x

if __name__ == '__main__':
    unittest.main()
