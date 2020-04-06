from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver
from gillespy2.core import log

def check_cpp_support():
    from gillespy2.solvers.cpp.example_models import Example
    try:
        model = Example()
        results = model.run(solver=SSACSolver)
        return True
    except Exception as e:
        log.warn('Unable to use C++ optimized SSA: {0}.  The performance of ' \
        'this package can be significantly increased if you install/configure GCC on ' \
        'this machine.'.format(e))
        return False

can_use_cpp = check_cpp_support()

__all__ = ['SSACSolver'] if can_use_cpp else []
