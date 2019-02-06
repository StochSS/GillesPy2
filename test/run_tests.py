import unittest, sys, os
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '..')))

import test_empty_model, test_simple_model

if __name__ == '__main__':
    modules = [
        test_empty_model,
        test_simple_model
        ]

    for module in modules:
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner()
        result = runner.run(suite)
    
    sys.exit(not result.wasSuccessful())