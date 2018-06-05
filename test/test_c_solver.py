import unittest
import tempfile
from gillespy2 import SSACSolver
from gillespyError import SolverError, DirectoryError, BuildError, ExecutionError
from gillespy2.example_models import Example

class TestSSACSolver(unittest.TestCase):
    def test_create(self):
        model = Example()
        solver = SSACSolver(model)
    
    def test_directory_exists(self):
        with self.assertRaises(DirectoryError):
            directory = tempfile.TemporaryDirectory()
            model = Example()
            solver = SSACSolver(model, directory)


if __name__ == '__main__':
    unittest.main()
