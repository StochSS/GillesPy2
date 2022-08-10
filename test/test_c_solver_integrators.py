import unittest
from gillespy2 import log
from gillespy2 import TauHybridCSolver
from gillespy2 import ODECSolver
from example_models import MichaelisMenten

class TestCSolverIntegrators(unittest.TestCase):
    """
    Test cases for ensuring that solvers with a dependency on SUNDIALS/CVODE integrators are configured correctly.
    """

    integrator_solvers = [
        TauHybridCSolver,
        ODECSolver,
    ]
    test_model = MichaelisMenten()

    @classmethod
    def setUpClass(cls):
        cls.integrator_solvers = [
            solver_class(model=cls.test_model) for solver_class in cls.integrator_solvers
        ]

    def test_tolerance_options_accepted(self):
        """
        Ensure that passing --rtol, --atol, etc. as integrator options is validated successfully.
        """
        integrator_options = {
            "rtol": 1e-4,
            "atol": 1e-4,
            "max_step": 1000,
        }
        for solver in self.integrator_solvers:
            with self.subTest(solver=solver.name), \
                    self.assertRaises(AssertionError, msg="Logger unexpectedly invoked when integrator_options passed"), \
                    self.assertLogs(log):
                solver.run(t=1, integrator_options=integrator_options)

if __name__ == "__main__":
    unittest.main(TestCSolverIntegrators)
