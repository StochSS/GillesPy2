'''
SpatialPy is a Python 3 package for simulation of
spatial deterministic/stochastic reaction-diffusion-advection problems
Copyright (C) 2021 SpatialPy developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU GENERAL PUBLIC LICENSE Version 3 as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU GENERAL PUBLIC LICENSE Version 3 for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import unittest
import os
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor


class TestNotebooks(unittest.TestCase):
    errors = {}
    ep = ExecutePreprocessor(timeout=600, kernel_name='python3', allow_errors=True)
    
    def test_notebooks(self):
        for root, dirs, files in os.walk("../examples/"):
            for file in files:
                if file.endswith(".ipynb"):
                     with open(os.path.join(root, file)) as f:
                        print('Executing {}...'.format(file))
                        nb = nbformat.read(f, as_version=nbformat.NO_CONVERT)
                        try:
                            ep.preprocess(nb, {'metadata': {'path': root}})
                        except Exception as err:
                            print('Error executing the notebook "{}".\n\n'.format(file))
                            errors[file] = err
        for fname, err in errors.items():
            if len(err.__str__()) > 500:
                print('{}:\n{}\n...\n{}'.format(fname, err.__str__()[:251], err.__str__()[-251:]))
            else:
                print('{}:\n{}'.format(fname, err))

