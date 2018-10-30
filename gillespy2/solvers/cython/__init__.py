from gillespy2.core import log
try:
    import pyximport
    import numpy as np
    pyximport.install(setup_args={'include_dirs': np.get_include()})
    from gillespy2.solvers.cython.cython_ssa_solver import CythonSSASolver
    can_use_cython = True
    log.debug("Successful Import of CythonSSASolver.")
except Exception as e:
    log.warn(" Unable to use Cython optimized SSA: {0}\nThe performance of this\
                 package can be significantly increased if you install Cython.".format(e))
    can_use_cython = False
