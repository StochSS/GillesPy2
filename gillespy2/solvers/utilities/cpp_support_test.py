
"""
This file contains a function and variable for testing a machines support of GillesPy2 C++ solvers.
Used in model.py
"""


def check_cpp_support():
    from gillespy2.solvers.cpp.example_models import Example
    from gillespy2 import SSACSolver
    try:
        model = Example()
        results = model.run(solver=SSACSolver, cpp_support=True)
        return True
    except Exception as e:
        log.warn('Unable to use C++ optimized SSA: {0}.  The performance of ' \
        'this package can be significantly increased if you install/configure GCC on ' \
        'this machine.'.format(e))
        return False


cpp_support = check_cpp_support()
