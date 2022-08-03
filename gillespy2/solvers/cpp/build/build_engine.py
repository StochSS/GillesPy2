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

import os
import shutil
import tempfile
import platform
from pathlib import Path
from typing import Union

import gillespy2
from gillespy2.core import Model, gillespyError

from . import template_gen
from .make import Make


class BuildEngine():
    template_definitions_name = "template_definitions.h"
    template_options_name = "template_opts.h"

    def __init__(self, debug: bool = False, output_dir: str = None):
        self.self_dir = Path(__file__).parent
        self.cpp_dir = self.self_dir.joinpath("../c_base").resolve()
        if platform.system() == 'Windows':
            self.makefile = self.cpp_dir.joinpath("Makefile.win32")
        else:
            self.makefile = self.cpp_dir.joinpath("Makefile")
        self.src_template_dir = self.cpp_dir.joinpath("template")
        self.output_dir = output_dir

        self.debug = debug

        # The output_dir will be generated if it does not exist on prepare().
        # Output files are all rooted relative to the output_dir.
        if self.output_dir is not None:
            self.output_dir = Path(output_dir)

            if not self.output_dir.is_dir():
                self.output_dir.mkdir(parents=True)

    def __get_cache_dir(self) -> Path:
        cache_dir = gillespy2._global_cache

        if not cache_dir.is_dir() or not os.access(str(cache_dir), os.R_OK):
            cache_dir = Path(os.path.expanduser("~"), ".gillespy2", "cache")

        return cache_dir

    @classmethod
    def get_missing_dependencies(cls) -> "list[str]":
        """
        Determine which dependencies are missing on the system, if any.

        :returns: A list of missing dependencies.
        """

        dependencies = ["g++", "make"]
        missing = [(dep) for dep in dependencies if shutil.which(dep) is None]

        return missing

    def prepare(self, model: "Union[Model, template_gen.SanitizedModel]", variable=False) -> str:
        """
        Prepare the template directory for compilation.
        The following operations will be performed:
            1. If no cached object files are found, prebuild them.
            2. Copy the C++ template directory into a new temp directory.
            3. Remove the sample template_definitions.h file.
            4. Generate and write a template_definitions.h file from the model.

        :param model: A GillesPy2 model who's template definitions will be generated.
        :type model: gillespy2.Model

        :param variable: A template_gen argument requirement which enables support for non-constant param values.
        :type variable: bool

        :returns: The path of the output directory.
        """

        # If the output directory is None, create and set it to a temporary directory.
        if self.output_dir is None:
            self.output_dir = Path(tempfile.mkdtemp(
                prefix='gillespy2_build_', dir=os.environ.get('GILLESPY2_TMPDIR'))
            )

        # If object files haven't been compiled yet, go ahead and compile them with make.
        # Precompilation only happens if the cache is enabled but hasn't been built yet.
        # Make target for individual simulation will still succeed if prebuild() isn't called.
        self.obj_dir = self.output_dir.joinpath("gillespy2_obj")
        self.template_dir = self.output_dir.joinpath("gillespy2_template")

        if gillespy2.cache_enabled:
            self.obj_dir = self.__get_cache_dir()

        if not self.obj_dir.is_dir():
            self.obj_dir.mkdir(parents=True)

        # Copy the C++ template directory to the temp directory.
        shutil.copytree(self.src_template_dir, self.template_dir)

        # If a raw GillesPy2 model was provided, convert it to a sanitized model.
        if isinstance(model, gillespy2.Model):
            model = template_gen.SanitizedModel(model, variable=variable)
        elif not isinstance(model, template_gen.SanitizedModel) and type(model).__name__ == "SanitizedModel":
            raise TypeError(f"Build engine expected gillespy2.Model or SanitizedModel type: received {type(model)} , __name__={type(model).__name__}")

        # Build the template and write it to the temp directory and remove the sample template_definitions header.
        template_file = self.template_dir.joinpath(self.template_definitions_name)
        template_file.unlink()
        template_gen.write_definitions(str(template_file), model.get_template())
        custom_definitions = model.get_options()
        if custom_definitions is not None:
            options_file = self.template_dir.joinpath(self.template_options_name)
            options_file.unlink()
            template_gen.write_definitions(str(options_file), custom_definitions)

        # With all required information gathered, create a Make instance.
        self.make = Make(str(self.makefile), str(self.output_dir), str(self.obj_dir))

        return self.output_dir

    def build_cache(self, cache_dir: str, force_rebuild: bool = False):
        """
        Build object dependencies and cache into directory for later use.

        :param cache_dir: The directory to use as a cache.
        :type cache_dir: str

        :param force_rebuild: Delete and rebuild the cache directory.
        :type force_rebuild: bool
        """

        make = Make(self.makefile, cache_dir, cache_dir)
        make.prebuild()

    def build_simulation(self, simulation_name: "str", definitions: "dict[str, str]" = None) -> str:
        """
        Build the solver to the temp directory.

        :param simulation_name: The name of the simulation to build. For example, ODESimulation.
        :type simulation_name: str

        :param definitions: Dictionary of environment variables to be overriden when the Makefile is invoked.
        Intended for use when running through a debugger environment or profiler.
        :type definitions: dict[str, str]

        :return: The path of the newly build solver executable.
        """
        if definitions is None:
            definitions = {}

        if self.make is None:
            raise gillespyError.BuildError(
                "Failed to build the simulation. The build environment has not yet been prepared.\n"
                "To fix, call `BuildEngine.prepare()` prior to attempting to build the simulation."
            )

        self.make.build_simulation(simulation_name, template_dir=str(self.template_dir), **definitions)
        return str(self.make.output_file)

    def get_executable_path(self) -> str:
        """
        Resolves the filepath of the simulation executable.
        Only valid after the simulation has been built.

        :return: String containing path to executable.
            None if no executable exists.
        """
        if not os.path.exists(self.make.output_file):
            return None
        return str(self.make.output_file)

    def clean(self):
        """
        Delete the output directory and all other associated build artifacts.
        """

        if self.debug:
            return

        if self.output_dir is None:
            return

        if self.output_dir.exists():
            shutil.rmtree(self.output_dir, ignore_errors=True)
