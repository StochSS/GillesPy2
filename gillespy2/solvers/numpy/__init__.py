from gillespy2.core import log
try:
    import numpy as np
    can_use_numpy = True
    from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
    from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver
    from gillespy2.solvers.numpy.basic_tau_leaping_solver import BasicTauLeapingSolver
    from gillespy2.solvers.numpy.basic_tau_hybrid_solver import BasicTauHybridSolver
    log.debug("Successful Import of NumPy solvers.")
except Exception as e:
    log.warn(" Unable to use NumPy: {0}. The performance of this package can be significantly increased if you install NumPy.".format(e))
    can_use_numpy = False


__all__ = ['NumPySSASolver', 'BasicODESolver', 'BasicTauLeapingSolver', 'BasicTauHybridSolver'] if can_use_numpy else []
