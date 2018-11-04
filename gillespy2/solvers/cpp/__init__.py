from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.solvers.cpp.test_ssa_c_solver import TestSSACSolver


def check_cpp_support():
    import unittest
    import os
    with open(os.devnull, 'w') as null_stream:
        runner = unittest.TextTestRunner(failfast=True, stream=null_stream)
        result = runner.run(unittest.makeSuite(TestSSACSolver))
    return result.testsRun > 0 and len(result.errors) == 0


can_use_cpp = check_cpp_support()

__all__ = ['SSACSolver'] if can_use_cpp else []
