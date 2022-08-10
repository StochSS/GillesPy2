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

from pathlib import Path
import unittest
import shutil
import os
from example_models import create_dimerization

import gillespy2
from gillespy2.solvers.cpp.build.build_engine import BuildEngine


class TestBuildEngine(unittest.TestCase):
    test_model = create_dimerization()

    # List of solver names.
    # Each entry represents the Makefile target for a C++ solver.
    solver_names = [
        "ssa",
        "ode",
        "tau_leap",
    ]

    def setUp(self) -> None:
        # Disable the gillespy2 cache.
        gillespy2.cache_enabled = False

        self.build_engine = BuildEngine()
        self.tmp_dir = self.build_engine.prepare(self.test_model, variable=False)

    def tearDown(self) -> None:
        self.build_engine.clean()

        if self.tmp_dir is not None and os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir, ignore_errors=True)

    def test_default_layout(self):
        """
        Ensure that with the default, cache-less output, the template and obj files
          are subdirectories of the temp directory.
        This is to ensure that they get cleaned up when the temp directory is cleaned up.
        """
        # This test should be attempted with each solver target.
        for solver_name in self.solver_names:
            with self.subTest(solver=solver_name):
                self.build_engine.build_simulation(solver_name)

        template_dir = self.build_engine.template_dir
        obj_dir = self.build_engine.obj_dir

        self.assertTrue(str(template_dir.resolve()).startswith(str(Path(self.tmp_dir).resolve())),
                        "Template directory was not placed in temp directory when cache is disabled")
        self.assertTrue(str(obj_dir.resolve()).startswith(str(Path(self.tmp_dir).resolve())),
                        "Dependency obj directory was not placed in temp directory when cache is disabled")

        # Ensure that the dependencies are cleaned up when .clean() is called.
        self.build_engine.clean()
        self.assertFalse(template_dir.exists(),
                         "Template directory was not cleaned up when cache is disabled")
        self.assertFalse(obj_dir.exists(),
                         "Dependency obj directory was not cleaned up when cache is disabled")

    def test_clean(self):
        """
        Ensure that with the default, cache-less output, the output files are cleaned up.
        """
        # Test clean-up with manual call to .clean()
        self.build_engine.clean()
        self.assertFalse(os.path.exists(self.tmp_dir), "Build engine output directory not cleaned up.")

    def test_template_file(self):
        """
        Ensure that the resulting template file exists and was prepared correctly.
        """
        template_file = self.build_engine.template_dir.joinpath(
                        self.build_engine.template_definitions_name)
        self.assertTrue(template_file.exists(), "Simulation template header file could not be found")

        # Test to ensure the contents of the output template_file contains ONLY macro definitions,
        # nothing else.
        with open(template_file) as template:
            for line in template.readlines():
                self.assertTrue(line.startswith("#"), "Output template file contains an invalid definition.")

    def test_build_output(self):
        """
        Ensure that the default, cache-less output has expected behavior.
        Should result in an expected definitions file and an executable simulation.
        """
        simulation_file = self.build_engine.output_dir.joinpath(
            "GillesPy2_Simulation.exe" if os.name == "nt" else "GillesPy2_Simulation.out")

        for solver_name in self.solver_names:
            with self.subTest(solver=solver_name):
                self.build_engine.build_simulation(solver_name)

        self.assertTrue(simulation_file.exists(), "Simulation executable could not be found")
        self.assertTrue(os.access(str(simulation_file), os.X_OK), "Simulation output is not executable")


if __name__ == '__main__':
    unittest.main()
