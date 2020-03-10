from gillespy2.core import log
try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs': np.get_include()})
    from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver
    can_use_cython = True
    log.debug("Successful Import of Cython solvers.")
except Exception as e:
    can_use_cython = False

__all__ = ['CythonSSASolver'] if can_use_cython else []
