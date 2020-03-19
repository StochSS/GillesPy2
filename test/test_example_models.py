import unittest
from gillespy2.example_models import *
from gillespy2.core.gillespyError import *
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver

class TestExampleModels(unittest.TestCase):

    def test_trichloroethylene_example(self):
        trichloroethylene_model = Trichloroethylene()
        results = trichloroethylene_model.run()

    def test_lacOperon_example(self):
        lacOperon_model = LacOperon()
        results = lacOperon_model.run(solver=BasicODESolver)

    def test_schlogl_example(self):
        schlogl_model = Schlogl()
        results = schlogl_model.run()

    def test_michaelisMenten_example(self):
        michaelisMenten_model = MichaelisMenten()
        results = michaelisMenten_model.run()

    def test_toggleSwitch_example(self):
        toggleSwitch_model = ToggleSwitch()
        results = toggleSwitch_model.run()

    def test_example_example(self):
        example_model = Example()
        results = example_model.run()

    def test_tyson2StateOscillator_example(self):
        tyson2StateOscillator_model = Tyson2StateOscillator()
        results = tyson2StateOscillator_model.run()


if __name__ == '__main__':
    unittest.main()
