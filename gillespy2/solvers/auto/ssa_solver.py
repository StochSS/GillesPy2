from gillespy2.core import log
from gillespy2.solvers.python import BasicSSASolver


def get_best_ssa_solver():
    ssa_solver = BasicSSASolver

    from gillespy2.solvers.numpy import can_use_numpy
    if can_use_numpy:
        from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
        log.debug("Successful Import of NumPySSASolver.")
        ssa_solver = NumPySSASolver

    from gillespy2.solvers.cython import can_use_cython
    if can_use_cython:
        from gillespy2.solvers.cython import CythonSSASolver
        log.debug("Successful Import of CythonSSASolver.")
        ssa_solver = CythonSSASolver
    return ssa_solver


SSASolver = get_best_ssa_solver()
