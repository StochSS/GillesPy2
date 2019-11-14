import unittest, sys, os
import argparse

try:
    import pyximport
    pyximport.install()
except Exception:
    pass

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mode', default='develop', choices=['develop', 'release'],
                    help='Run tests in develop mode or release mode.')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.mode == 'develop':
        print('Running tests in develop mode. Appending repository directory to system path.')
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


    import test_cython_ssa_solver
    import test_empty_model
    import test_model
    import test_ode_solver
    import test_tau_hybrid_solver
    import test_simple_model
    import test_ssa_c_solver
    import test_SBML
    import test_all_solvers

    modules = [
        test_cython_ssa_solver,
        test_empty_model,
        test_model,
        test_ode_solver,
        test_tau_hybrid_solver,
        test_simple_model,
        test_ssa_c_solver,
        test_SBML,
        test_all_solvers
    ]

    for module in modules:
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner(failfast=args.mode == 'develop')

        print("Executing: {}".format(module))
        result = runner.run(suite)
        print('=' * 70)
        if not result.wasSuccessful():
            sys.exit(not result.wasSuccessful())
