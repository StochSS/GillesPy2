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

"""
Test to determine if the behavior exhibited by the new check_cpp_support() function matches
what is expected by the previous try/catch implementation.
"""

import unittest
from gillespy2.solvers.cpp.build.build_engine import BuildEngine
from example_models import create_decay
from gillespy2 import SSACSolver

class TestCheckCPPSupport(unittest.TestCase):
    def test_check_cpp_support(self):
        self.assertEqual(not len(BuildEngine.get_missing_dependencies()), self.old_check_cpp_support())

    def old_check_cpp_support(self):
        try:
            model = create_decay()
            solver = SSACSolver(model=model)
            results = model.run(solver=solver)
            return True
        except Exception as e:
            from gillespy2.core import log
            log.warn('Unable to use C++ optimized SSA: {0}.  The performance of ' \
            'this package can be significantly increased if you install/configure GCC on ' \
            'this machine.'.format(e))
            return False

if __name__ == '__main__':
    unittest.main()
