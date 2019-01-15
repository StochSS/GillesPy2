import unittest

loader = unittest.TestLoader()
start_dir = '.'
suite = loader.discover(start_dir)

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite)