# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

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
import os
import sys
import unittest
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
        print('Running unit tests in develop mode. Appending repository directory to system path.')
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    from unit_tests import test_species
    from unit_tests import test_parameter
    from unit_tests import test_reaction
    from unit_tests import test_raterule
    from unit_tests import test_assignmentrule
    from unit_tests import test_event_assignment
    from unit_tests import test_timespan
    from unit_tests import test_model

    catagories = {
        "State Components": [test_species, test_parameter, test_timespan],
        "Action Dependents": [test_event_assignment],
        "Action Components": [test_reaction, test_raterule, test_assignmentrule],
        "Model": [test_model]
    }

    for name, modules in catagories.items():
        print(f"Running unit tests for {name}")
        for module in modules:
            suite = unittest.TestLoader().loadTestsFromModule(module)
            runner = unittest.TextTestRunner(failfast=args.mode == 'develop')

            print(f"Executing: {module}")
            result = runner.run(suite)
            print('=' * 70)
            if not result.wasSuccessful():
                sys.exit(not result.wasSuccessful())
