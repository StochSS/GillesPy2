import shutil, tempfile
from pathlib import Path

from . import template_gen
from .make import Make
from gillespy2.core import Model

class BuildEngine():
    template_definitions_name = "template_definitions.h"

    def __init__(self, debug: bool = False, temp_dir: str = None):
        self.self_dir = Path(__file__).parent
        self.cpp_dir = self.self_dir.joinpath("../c_base").resolve()
        self.template_dir = self.cpp_dir.joinpath("template")
        self.makefile = self.cpp_dir.joinpath("Makefile")
        self.temp_dir = Path(temp_dir)

        self.debug = debug

    def __enter__(self):
        # If the temp dir is not set, create one.
        if self.temp_dir is None:
            self.temp_dir = Path(tempfile.mktemp())

        # If the temp dir is set ensure that it exists, if not, create it.
        elif not self.temp_dir.is_dir():
            self.temp_dir.mkdir()

        self.make = Make(self.makefile, self.temp_dir)
        return self

    def __exit__(self, type, value, traceback):
        """
        Delete the temp directory and all other associated build artifacts.
        """

        # If debug is True do nothing.
        if self.debug:
            return

        shutil.rmtree(self.temp_dir)

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

        # Copy the C++ template directory to the temp directory and remove the sample template_definitions header.
        shutil.copytree(self.template_dir, self.temp_dir, dirs_exist_ok=True)
        self.temp_dir.joinpath(self.template_definitions_name).unlink()

        # Build the template and write it to the temp directory.
        template_gen.write_template(self.temp_dir.joinpath(self.template_definitions_name), model, variable)

    def build_simulation(self, simulation_name: str) -> str:
        """
        Build the solver to the temp directory.

        :param simulation_name: The name of the simulation to build. For example, ODESimulation.
        :type simulation_name: str

        :return: The path of the newly build solver executable.
        """

        self.make.build_solver(simulation_name, template_dir=self.template_dir)