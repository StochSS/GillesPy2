import unittest
from .example_models import MichaelisMenten
from gillespy2.solvers.cpp import SSACSolver, ODECSolver
from .perf.gprof import run_profiler

class MyTestCase(unittest.TestCase):
    def test_profiler(self):
        model = self.test_models[0]
        solver = self.test_solvers[0]

        time = run_profiler(model, solver(model=model))
        print(f"Time: {time}")

    def setUp(self) -> None:
        self.test_models = [
            MichaelisMenten(),
        ]
        self.test_solvers = [
            SSACSolver,
        ]

if __name__ == '__main__':
    unittest.main()
