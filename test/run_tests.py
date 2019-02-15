import unittest, sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import test_basic_tau_hybrid_solver
import test_basic_tau_leaping_solver
# import test_cython_ssa_solver
import test_empty_model
import test_model
import test_ode_solver
import test_simple_model
import test_ssa_solver
import test_ssa_c_solver

if __name__ == '__main__':
    modules = [
        test_basic_tau_hybrid_solver,
        test_basic_tau_leaping_solver,
        # test_cython_ssa_solver,
        test_empty_model,
        test_model,
        test_ode_solver,
        test_simple_model,
        test_ssa_solver,
        test_ssa_c_solver
    ]

    for module in modules:
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner()

        print("Executing: {}".format(module))
        result = runner.run(suite)
        print('=' * 70)
        if not result.wasSuccessful():
            sys.exit(not result.wasSuccessful())
