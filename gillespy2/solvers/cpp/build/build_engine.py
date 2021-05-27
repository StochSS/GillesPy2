import shutil, tempfile
import os
from pathlib import Path

from pexpect import which

from . import template_gen
from .make import Make
from gillespy2.core import Model

class BuildEngine():
    template_definitions_name = "template_definitions.h"

    def __init__(self, debug: bool = False, output_dir: str = None):
        # If the temp_dir is None, make one. Else, ensure it exists.
        # Output files are all rooted relative to the temp_dir.
        if output_dir is None:
            self.temp_dir = Path(tempfile.mkdtemp())
        else:
            self.temp_dir = Path(output_dir)

        self.self_dir = Path(__file__).parent
        self.cpp_dir = self.self_dir.joinpath("../c_base").resolve()
        self.template_dir = self.temp_dir.joinpath("template")
        self.makefile = self.cpp_dir.joinpath("Makefile")

        self.debug = debug

        if not self.temp_dir.is_dir():
            self.temp_dir.mkdir()

        self.make = Make(str(self.makefile), str(self.temp_dir))

    @classmethod
    def get_missing_dependencies(cls):
        """
        Determine which dependencies are missing on the system, if any.

        :returns: A list of missing dependencies.
        """

        dependencies = ["g++", "make"]
        missing = [(dep) for dep in dependencies if shutil.which(dep) is None]

        return missing

    def prepare(self, model: Model, variable=False):
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
        """

        # If object files haven't been compiled yet, go ahead and compile them with make.
        self.make.prebuild()

        # Output files are rooted relative to the output directory.
        # Template files get placed in <output>/template.
        src_template_dir = self.cpp_dir.joinpath("template")

        # Copy the C++ template directory to the temp directory.
        shutil.copytree(src_template_dir, self.template_dir)

        # Build the template and write it to the temp directory and remove the sample template_definitions header.
        template_file = self.template_dir.joinpath(self.template_definitions_name)
        template_file.unlink()
        template_gen.write_template(str(template_file), model, variable)

    def build_simulation(self, simulation_name: str) -> str:
        """
        Build the solver to the temp directory.

        :param simulation_name: The name of the simulation to build. For example, ODESimulation.
        :type simulation_name: str

        :return: The path of the newly build solver executable.
        """

        self.make.build_simulation(simulation_name, template_dir=str(self.template_dir))
        return str(self.make.output_file)

    def clean(self):
        """
        Delete the output directory and all other associated build artifacts.
        """

        if self.debug:
            return

        shutil.rmtree(self.temp_dir)