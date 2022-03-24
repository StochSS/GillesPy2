
"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import unittest
import os
import tempfile
import sys
import copy
import numpy
import shutil


class TestCompileWSpaces(unittest.TestCase):
    
    def test_compile_w_spaces(self):
        prefix_base_dir = tempfile.mkdtemp()
        os.mkdir(prefix_base_dir+'/A SPACE')
        shutil.copytree(os.path.abspath(os.path.dirname(__file__))+"/../gillespy2",prefix_base_dir+'/A SPACE/gillespy2')
        old_path = copy.copy(sys.path)
        sys.path.insert(0,prefix_base_dir+'/A SPACE/')
        import gillespy2
        # 
        try:
            # create a model
            model = gillespy2.Model(name="test_compile_model")
            model.add_species([
                gillespy2.Species(name='A', initial_value=0),
                gillespy2.Species(name='B', initial_value=0)
            ])
            model.add_parameter([
                gillespy2.Parameter(name='k1', expression=1),
                gillespy2.Parameter(name='k2', expression=10)
            ])
            model.add_rate_rule([
                gillespy2.RateRule(name='Brate', variable='B', formula="cos(t)")
            ])
            model.add_reaction([
                gillespy2.Reaction(reactants={'A': 1}, products={}, 
                                   propensity_function="k1*B"),
                gillespy2.Reaction(reactants={}, products={'B': 1}, 
                                   rate='k2')
            ])
            model.timespan(numpy.array([0., 5., 10.]))                

            # run the model
            result = model.run()
        finally:
            # Cleanup
            sys.path = old_path
            shutil.rmtree(prefix_base_dir)
        # Checks
        self.assertTrue(gillespy2.__file__ == prefix_base_dir+'/A SPACE/gillespy2/__init__.py')


if __name__ == '__main__':
    unittest.main()

