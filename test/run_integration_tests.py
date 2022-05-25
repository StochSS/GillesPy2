# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest, sys, os
import argparse

try:
    import pyximport
    pyximport.install()
except Exception:
    pass

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mode', default='develop', choices=['develop', 'release'],
                    help='Run tests in develop mode or release mode.')

if __name__ == '__main__':
    args = parser.parse_args()
    if args.mode == 'develop':
        print('Running tests in develop mode. Appending repository directory to system path.')
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    import test_empty_model
    import test_build_engine
    import test_c_solvers
    import test_model
    import test_ode_solver
    import test_ode_c_solver
    import test_hybrid_solver
    import test_simple_model
    import test_ssa_c_solver
    import test_variable_solvers
    import test_tau_leaping_c_solver
    import test_tau_hybrid_c_solver
    import test_SBML
    import test_StochML
    import test_example_models
    import test_all_solvers
    import test_sys_init
    import test_results
    import test_propensity_parser
    import test_pause_resume
    import test_check_cpp_support
    import test_jsonify
    import test_notebooks

    modules = [
        test_empty_model,
        test_build_engine,
        test_c_solvers,
        test_model,
        test_ode_solver,
        test_ode_c_solver,
        test_tau_leaping_c_solver,
        test_tau_hybrid_c_solver,
        test_hybrid_solver,
        test_simple_model,
        test_ssa_c_solver,
        test_variable_solvers,
        test_pause_resume,
        test_SBML,
        test_StochML,
        test_example_models,
        test_all_solvers,
        test_sys_init,
        test_results,
        test_propensity_parser,
        test_check_cpp_support,
        test_jsonify,
        test_notebooks
    ]

    for module in modules:
        suite = unittest.TestLoader().loadTestsFromModule(module)
        runner = unittest.TextTestRunner(failfast=args.mode == 'develop')

        print("Executing: {}".format(module))
        result = runner.run(suite)
        print('=' * 70)
        if not result.wasSuccessful():
            sys.exit(not result.wasSuccessful())
