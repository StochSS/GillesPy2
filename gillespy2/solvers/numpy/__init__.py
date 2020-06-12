import logging
_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
_handler = logging.StreamHandler()
_handler.setFormatter(_formatter)
log = logging.getLogger()
log.setLevel(logging.WARN)
log.addHandler(_handler)

try:
    import numpy as np
    can_use_numpy = True
    from .ssa_solver import NumPySSASolver
    from .ode_solver import ODESolver
    from .tau_leaping_solver import TauLeapingSolver
    from .tau_hybrid_solver import TauHybridSolver
    log.debug("Successful Import of NumPy solvers.")

except Exception as e:
    log.warn(" Unable to use NumPy: {0}. The performance of this package can be significantly increased if you install NumPy.".format(e))
    can_use_numpy = False


__all__ = ['NumPySSASolver', 'ODESolver', 'TauLeapingSolver', 'TauHybridSolver'] if can_use_numpy else []
