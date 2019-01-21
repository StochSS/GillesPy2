import unittest
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))

loader = unittest.TestLoader()
start_dir = '.'
suite = loader.discover(start_dir)

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite)

