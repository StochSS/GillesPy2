"""
This file contains a function and variable for testing a machines support of GillesPy2 C++ solvers.
Used in model.py
"""

from gillespy2.solvers.cpp.build.build_engine import BuildEngine

def check_cpp_support():
    return not len(BuildEngine.get_missing_dependencies()) > 0

cpp_support = check_cpp_support()
