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
''' Unit tests module for gillespy2.Species. '''
import time
import unittest
from datetime import datetime

from gillespy2 import Species
from gillespy2 import SpeciesError

class TestSpecies(unittest.TestCase): # pylint: disable=too-many-public-methods
    ''' Unit tests class for gillespy2.Species. '''
    def setUp(self):
        self.species = Species(name="test_species", initial_value=5)

    def test_constructor(self):
        """ Test the Species constructor. """
        species = Species(name="test_species", initial_value=5)
        self.assertEqual(species.name, "test_species")
        self.assertEqual(species.initial_value, 5)

    def test_constructor__no_name(self):
        """ Test the Species constructor without name. """
        with self.assertRaises(SpeciesError):
            Species(initial_value=0)

    def test_constructor__invalid_name(self):
        """ Test the Species constructor with non-str or empty str name. """
        test_names = ["", None, 0]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(SpeciesError):
                    Species(name=test_name, initial_value=0)

    def test_constructor__no_initial_value(self):
        """ Test the Species constructor without initial_value. """
        species = Species(name="test_species")
        self.assertEqual(species.initial_value, 0)

    def test_constructor__negative_initial_value(self):
        """ Test the Species constructor with negative initial_value. """
        with self.assertRaises(SpeciesError):
            Species(name="test_species", initial_value=-1)

    def test_constructor__negative_initial_value_allowed(self):
        """
        Test the Species constructor with negative initial_value and allow_negative_populations.
        """
        species = Species(
            name="test_species", initial_value=-1, allow_negative_populations=True
        )
        self.assertEqual(species.initial_value, -1)

    def test_constructor__invalid_initial_value(self):
        """ Test the Species constructor with non-int or non-float initial_value. """
        test_ivs = [None, "0%"]
        for test_iv in test_ivs:
            with self.subTest(initial_value=test_iv):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", initial_value=test_iv)

    def test_constructor__invalid_constant(self):
        """ Test the Species constructor with non-bool constant. """
        test_constants = [None, "True"]
        for test_constant in test_constants:
            with self.subTest(constant=test_constant):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", constant=test_constant)

    def test_constructor__invalid_boundary_condition(self):
        """ Test the Species constructor with non-bool boundary_condition. """
        test_bcs = [None, "True"]
        for test_bc in test_bcs:
            with self.subTest(boundary_condition=test_bc):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", boundary_condition=test_bc)

    def test_constructor__invalid_mode(self):
        """ Test the Species constructor with and invalid mode. """
        with self.assertRaises(SpeciesError):
            Species(name="test_species", mode="deterministic")

    def test_constructor__continuous_mode_float_initial_value(self):
        """ Test the Species constructor with continuous mode and float initial value. """
        species = Species(name="test_species", initial_value=0.1, mode="continuous")
        self.assertEqual(species.mode, "continuous")
        self.assertIsInstance(species.initial_value, float)
        self.assertEqual(species.initial_value, 0.1)

    def test_constructor__continuous_mode_int_initial_value(self):
        """ Test the Species constructor with continuous mode and int initial value. """
        species = Species(name="test_species", initial_value=1, mode="continuous")
        self.assertEqual(species.mode, "continuous")
        self.assertIsInstance(species.initial_value, float)
        self.assertEqual(species.initial_value, 1.0)

    def test_constructor__non_continuous_mode_float_initial_value(self):
        """ Test the Species constructor with non-continuous mode and float initial value. """
        test_modes = ['discrete', 'dynamic', None]
        for test_mode in test_modes:
            with self.subTest(mode=test_mode):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", initial_value=0.1, mode=test_mode)

    def test_constructor__non_continuous_mode_int_initial_value(self):
        """ Test the Species constructor with non-continuous mode and int initial value. """
        test_modes = ['discrete', 'dynamic', None]
        for test_mode in test_modes:
            with self.subTest(mode=test_mode):
                species = Species(name="test_species", initial_value=1, mode=test_mode)
                self.assertEqual(species.mode, test_mode)
                self.assertIsInstance(species.initial_value, int)
                self.assertEqual(species.initial_value, 1)

    def test_constructor__invalid_allow_negative_populations(self):
        """ Test the Species constructor with non-bool allow_negative_populations. """
        test_anps = [None, "True"]
        for test_anp in test_anps:
            with self.subTest(allow_negative_populations=test_anp):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", allow_negative_populations=test_anp)

    def test_constructor__int_switch_tol(self):
        """ Test the Species constructor with int switch_tol. """
        species = Species(name="test_species", switch_tol=1)
        self.assertEqual(species.switch_tol, 1.0)

    def test_constructor__invalid_switch_tol(self):
        """ Test the Species constructor with non-int, -float switch_tol. """
        test_sts = [None, "1%", -1]
        for test_st in test_sts:
            with self.subTest(switch_tol=test_st):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", switch_tol=test_st)

    def test_constructor__int_switch_min(self):
        """ Test the Species constructor with int switch_min. """
        species = Species(name="test_species", switch_min=1)
        self.assertEqual(species.switch_min, 1.0)

    def test_constructor__invalid_switch_min(self):
        """ Test the Species constructor with non-int, -float switch_min. """
        test_sms = [None, "1%", -1]
        for test_sm in test_sms:
            with self.subTest(switch_min=test_sm):
                with self.assertRaises(SpeciesError):
                    Species(name="test_species", switch_min=test_sm)

    def test___str___(self):
        """ Test Species.__str__ method. """
        self.assertIsInstance(str(self.species), str)

    def test_set_initial_value__negative_initial_value(self):
        """ Test Species.set_initial_value with negative initial_value. """
        with self.assertRaises(SpeciesError):
            self.species.set_initial_value(-1)

    def test_set_initial_value__negative_initial_value_allowed(self):
        """
        Test Species.set_initial_value with negative initial_value and allow_negative_populations.
        """
        self.species.allow_negative_populations = True
        self.species.set_initial_value(-1)
        self.assertEqual(self.species.initial_value, -1)

    def test_set_initial_value__invalid_initial_value(self):
        """ Test Species.set_initial_value with non-int or non-float initial_value. """
        test_ivs = [None, "0%"]
        for test_iv in test_ivs:
            with self.subTest(initial_value=test_iv):
                with self.assertRaises(SpeciesError):
                    self.species.set_initial_value(test_iv)

    def test_set_initial_value__continuous_mode_float_initial_value(self):
        """ Test Species.set_initial_value with continuous mode and float initial value. """
        self.species.mode = "continuous"
        self.species.set_initial_value(0.1)
        self.assertIsInstance(self.species.initial_value, float)
        self.assertEqual(self.species.initial_value, 0.1)

    def test_set_initial_value__continuous_mode_int_initial_value(self):
        """ Test Species.set_initial_value with continuous mode and int initial value. """
        self.species.mode = "continuous"
        self.species.set_initial_value(1)
        self.assertIsInstance(self.species.initial_value, float)
        self.assertEqual(self.species.initial_value, 1.0)

    def test_set_initial_value__non_continuous_mode_float_initial_value(self):
        """ Test Species.set_initial_value with non-continuous mode and float initial value. """
        test_modes = ['discrete', 'dynamic', None]
        for test_mode in test_modes:
            with self.subTest(mode=test_mode):
                self.species.mode = test_mode
                with self.assertRaises(SpeciesError):
                    self.species.set_initial_value(0.1)

    def test_set_initial_value__non_continuous_mode_int_initial_value(self):
        """ Test Species.set_initial_value with non-continuous mode and int initial value. """
        test_modes = ['discrete', 'dynamic', None]
        for test_mode in test_modes:
            with self.subTest(mode=test_mode):
                self.species.mode = test_mode
                self.species.set_initial_value(1)
                self.assertIsInstance(self.species.initial_value, int)
                self.assertEqual(self.species.initial_value, 1)

    def test_validate__invalid_name(self):
        """ Test Species.validate with non-str or empty str name. """
        test_names = ["", None, 0]
        for test_name in test_names:
            with self.subTest(name=test_name):
                with self.assertRaises(SpeciesError):
                    self.species.name = test_name
                    self.species.validate()

    def test_validate__negative_initial_value(self):
        """ Test Species.validate with negative initial_value. """
        with self.assertRaises(SpeciesError):
            self.species.initial_value = -1
            self.species.validate()

    def test_validate__invalid_initial_value(self):
        """ Test Species.validate with non-int or non-float initial_value. """
        test_ivs = [None, "0%"]
        for test_iv in test_ivs:
            with self.subTest(initial_value=test_iv):
                with self.assertRaises(SpeciesError):
                    self.species.initial_value = test_iv
                    self.species.validate()

    def test_validate__invalid_constant(self):
        """ Test Species.validate with non-bool constant. """
        test_constants = [None, "True"]
        for test_constant in test_constants:
            with self.subTest(constant=test_constant):
                with self.assertRaises(SpeciesError):
                    self.species.constant = test_constant
                    self.species.validate()

    def test_validate__invalid_boundary_condition(self):
        """ Test Species.validate with non-bool boundary_condition. """
        test_bcs = [None, "True"]
        for test_bc in test_bcs:
            with self.subTest(boundary_condition=test_bc):
                with self.assertRaises(SpeciesError):
                    self.species.boundary_condition = test_bc
                    self.species.validate()

    def test_validate__invalid_mode(self):
        """ Test Species.validate with and invalid mode. """
        with self.assertRaises(SpeciesError):
            self.species.mode = "deterministic"
            self.species.validate()

    def test_validate__non_continuous_mode_float_initial_value(self):
        """ Test Species.validate with non-continuous mode and float initial value. """
        test_modes = ['discrete', 'dynamic', None]
        for test_mode in test_modes:
            with self.subTest(mode=test_mode):
                with self.assertRaises(SpeciesError):
                    self.species.initial_value = 0.5
                    self.species.mode = test_mode
                    self.species.validate()

    def test_validate__invalid_allow_negative_populations(self):
        """ Test Species.validate with non-bool allow_negative_populations. """
        test_anps = [None, "True"]
        for test_anp in test_anps:
            with self.subTest(allow_negative_populations=test_anp):
                with self.assertRaises(SpeciesError):
                    self.species.allow_negative_populations = test_anp
                    self.species.validate()

    def test_validate__invalid_switch_tol(self):
        """ Test Species.validate with non-int, -float switch_tol. """
        test_sts = [None, "1%", -1]
        for test_st in test_sts:
            with self.subTest(switch_tol=test_st):
                with self.assertRaises(SpeciesError):
                    self.species.switch_tol = test_st
                    self.species.validate()

    def test_validate__invalid_switch_min(self):
        """ Test Species.validate with non-int, -float switch_min. """
        test_sms = [None, "1%", -1]
        for test_sm in test_sms:
            with self.subTest(switch_min=test_sm):
                with self.assertRaises(SpeciesError):
                    self.species.switch_min = test_sm
                    self.species.validate()

    def test_comp_time_of_validate(self):
        """ Check the computation time of validate. """
        start = time.time()
        self.species.validate()
        tic = datetime.utcfromtimestamp(time.time() - start)
        msg = f"Total time to run validate on a species: {tic.strftime('%M mins %S secs %f msecs')}"
        print(f"\n<{'-'*88}>\n | {msg.ljust(84)} | \n<{'-'*88}>")
