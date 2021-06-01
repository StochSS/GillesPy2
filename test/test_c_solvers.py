import unittest
import numpy
import io
import os
import example_models
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder
from gillespy2.solvers.cpp import SSACSolver, ODECSolver, TauLeapingCSolver

class TestCSolvers(unittest.TestCase):
    """
    Run tests to ensure the setup process is valid for each C++ solver.
    Does not actually run any simulations, just validates their construction.
    Effectively an "integration test" for each C++ solver, making sure they
      play nice with the build engine and template gen.
    """

    test_model = example_models.Dimerization()
    solvers = {
        SSACSolver.target: SSACSolver(model=test_model),
        ODECSolver.target: ODECSolver(model=test_model),
        TauLeapingCSolver.target: TauLeapingCSolver(model=test_model),
    }
    solvers_variable = {
        SSACSolver.target: SSACSolver(model=test_model, variable=True),
        ODECSolver.target: ODECSolver(model=test_model, variable=True),
        TauLeapingCSolver.target: TauLeapingCSolver(model=test_model, variable=True),
    }
    
    def test_c_decoder(self):
        """
        Ensures that, given a certain input from stdout, the results will be
          properly formatted by the sim decoder.
        """
        example_input = """
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,3
        """.strip()
        mock_stdout = io.BytesIO(bytes(example_input, "utf-8"))
        mock_stdout = io.BufferedReader(mock_stdout)

        trajectories = numpy.zeros((2, 4, 4))
        reader = BasicSimDecoder(trajectories)
        reader.read(mock_stdout)
        result, time_stopped = reader.get_output()
        
        expected_result = numpy.array([
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ],
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ]
        ])

        self.assertEqual(time_stopped, 3)
        # Test both the passed-in trajectory list and the returned trajectories.
        self.assertTrue(numpy.all(result == expected_result))

    def test_solver_build(self):
        """
        Build each solver and ensure that they build properly.
        """
        # Test builds for non-variable solvers
        for solver_name, solver in self.solvers.items():
            with self.subTest(solver=solver_name):
                exe = solver._build(model=self.test_model,
                                    simulation_name=solver_name,
                                    variable=False,
                                    debug=False)

                self.assertTrue(os.path.isfile(exe),
                                "Built simulation output could not be found or is a directory.")
                self.assertTrue(os.access(exe, os.X_OK),
                                "Built simulation binaries are not executable.")
                solver.build_engine.clean()
                self.assertFalse(os.path.exists(exe),
                                 "Simulation output not cleaned up after call to .clean().")

        # Test builds for variable solvers
        for solver_name, solver in self.solvers_variable.items():
            with self.subTest(solver=solver_name):
                exe = solver._build(model=self.test_model,
                                    simulation_name=solver_name,
                                    variable=True,
                                    debug=False)
                self.assertTrue(os.path.isfile(exe),
                                "Built simulation output could not be found or is a directory.")
                self.assertTrue(os.access(exe, os.X_OK),
                                "Built simulation binaries are not executable.")
                solver.build_engine.clean()
                self.assertFalse(os.path.exists(exe),
                                 "Simulation output not cleaned up after call to .clean().")
