import os, tempfile, shutil

from . import make, template_gen
from gillespy2.core import Model

class BuildEngine():
    template_definitions_name = "template_definitions.h"

    def __init__(self):
        self.self_dir = os.path.dirname(os.path.abspath(__file__))
        self.cpp_dir = os.path.join(self.self_dir, "../c_base/")

    def __enter__(self):
        """
        Create a new temp directory.
        """

        self.temp_dir = tempfile.gettempdir()
        os.mkdir(self.temp_dir)

        return self

    def __exit__(self):
        """
        Delete the temp directory and all other associated build artifacts.
        """

        shutil.rmtree(self.temp_dir)

    def prepare(self):
        """
        Copy everything but the template definitions header to the temp directory.
        """

        shutil.copytree(self.cpp_dir, self.temp_dir)
        os.remove(os.path.join(self.temp_dir, self.template_definitions_name))

    def write_template(self, model: Model, variable=False):
        """
        Write model definitions to the template file. If the file does not exist, one is created.

        :param model: The model who's definitions will be written.
        :type model: gillespy2.Model

        :param variable: Set true to allow for non-constant parameter values.
        :type variable: bool
        """

        template_gen.write_template(os.path.join(self.temp_dir, self.template_definitions_name), model, variable)