from gillespy2.core import log
from gillespy2.solvers.python import BasicSSASolver


ssa_solvers = []


def get_best_ssa_solver():
    global ssa_solvers
    ssa_solvers = [BasicSSASolver]

    from gillespy2.solvers.numpy import can_use_numpy
    if can_use_numpy:
        from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
        log.debug("Successful Import of NumPySSASolver.")
        ssa_solvers.append(NumPySSASolver)

    from gillespy2.solvers.cython import can_use_cython
    if can_use_cython:
        from gillespy2.solvers.cython import CythonSSASolver
        log.debug("Successful Import of CythonSSASolver.")
        ssa_solvers.append(CythonSSASolver)

    import unittest
    from gillespy2.solvers.cpp.test_ssa_c_solver import TestSSACSolver
    import os
    with open(os.devnull, 'w') as null_stream:
        runner = unittest.TextTestRunner(failfast=True, stream=null_stream)
        result = runner.run(unittest.makeSuite(TestSSACSolver))
    if result.testsRun > 0 and len(result.errors) == 0:
        from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
        ssa_solvers.append(SSACSolver)

    return ssa_solvers[-1]


SSASolver = get_best_ssa_solver()
