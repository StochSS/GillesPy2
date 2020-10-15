
"""
This file contains a function and variable for testing a machines support of GillesPy2 C++ solvers.
Used in model.py
"""

def check_cpp_support():
    import shutil

    dependencies = ['g++', 'make']
    missing = []
    any_missing = False

    for dependency in dependencies:
        if shutil.which(dependency) != None:
            continue

        missing.append(dependency)
        any_missing = True

    if any_missing is True:
        from gillespy2.core import log
        log.warn('Unable to use C++ optimized SSA due to one or more missing dependencies: {0}. '
        'The performance of this package can be significantly increased if you install/configure '
        'these on your machine.'.format(missing))

    return not any_missing

cpp_support = check_cpp_support()
