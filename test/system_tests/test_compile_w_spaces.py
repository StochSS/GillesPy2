
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

import os
import sys
import copy
import numpy
import shutil
import unittest
import tempfile
import platform



class TestCompileWSpaces(unittest.TestCase):
    
    def setUp(self):
        self.prefix_base_dir = tempfile.mkdtemp()
        os.mkdir(self.prefix_base_dir+'/A SPACE')
        shutil.copytree(os.path.abspath(os.path.dirname(__file__)+"/../../gillespy2"),self.prefix_base_dir+'/A SPACE/gillespy2')
        self.old_path = copy.copy(sys.path)
        sys.path.insert(0,self.prefix_base_dir+'/A SPACE/')
        
    def tearDown(self):
        sys.path = self.old_path
        shutil.rmtree(self.prefix_base_dir)

    def test_compile_w_spaces(self):
        import gillespy2
        self.solvers = [
            gillespy2.ODECSolver,
            gillespy2.SSACSolver,
            gillespy2.TauLeapingCSolver,
            gillespy2.TauHybridCSolver
        ]

        # create a model
        model = Model(name="Michaelis_Menten")

        # parameters
        rate1 = Parameter(name='rate1', expression=0.0017)
        rate2 = Parameter(name='rate2', expression=0.5)
        rate3 = Parameter(name='rate3', expression=0.1)
        model.add_parameter([rate1, rate2, rate3])

        # Species
        A = Species(name='A', initial_value=301)
        B = Species(name='B', initial_value=120)
        C = Species(name='C', initial_value=0)
        D = Species(name='D', initial_value=0)
        model.add_species([A, B, C, D])

        # reactions
        r1 = Reaction(name="r1", reactants={A: 1, B: 1}, products={C: 1}, rate=rate1)

        r2 = Reaction(name="r2", reactants={C: 1}, products={A: 1, B: 1}, rate=rate2)

        r3 = Reaction(name="r3", reactants={C: 1}, products={B: 1, D: 1}, rate=rate3)
        model.add_reaction([r1, r2, r3])
        model.timespan(np.linspace(0, 100, 101))

        # run the model
        for solver in self.solvers:
            with self.subTest(solver=solver.name):
                if platform.system() == "Windows":
                    with self.assertRaises(gillespy2.SimulationError):
                        model.run(solver=solver)
                else:
                    result = model.run(solver=solver)
                    self.assertTrue(gillespy2.__file__ == os.path.join(self.prefix_base_dir, 'A SPACE/gillespy2/__init__.py'))

if __name__ == '__main__':
    unittest.main()
