from gillespy2.core import log

def get_best_ssa_solver(omit_cpp=False, omit_cython=True, omit_numpy=False):
    if not omit_cpp:
        from gillespy2.solvers.cpp import can_use_cpp
        if can_use_cpp:
            from gillespy2.solvers.cpp import SSACSolver
            log.debug("Successful Import of SSACSolver.")
            return SSACSolver

    if not omit_cython:
        from gillespy2.solvers.cython import can_use_cython
        if can_use_cython:
            from gillespy2.solvers.cython import CythonSSASolver
            log.debug("Successful Import of CythonSSASolver.")
            return CythonSSASolver

    if not omit_numpy:
        from gillespy2.solvers.numpy import can_use_numpy
        if can_use_numpy:
            from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
            log.debug("Successful Import of NumPySSASolver.")
            return NumPySSASolver

    log.warn('Minimum software requirements not met.  Please install Numpy.')
    return None


SSASolver = get_best_ssa_solver()


__all__ = ['SSASolver', 'get_best_ssa_solver']
