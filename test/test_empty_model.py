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
from gillespy2.core.model import Model


class EmptyModel(Model):
    def __init__(self, parameter_values=None):
        Model.__init__(self)


class TestEmptyModel(unittest.TestCase):
    def setUp(self):
        self.model = EmptyModel()

    def test_model_creation(self):
        name = self.model.name
        self.assertEqual(name, "", msg="Unexpected value: {}".format(name))


if __name__ == '__main__':
    unittest.main()
