import unittest
from unittest.case import expectedFailure
import numpy
import io
import example_models
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder
from gillespy2.solvers.cpp import SSACSolver, ODECSolver, TauLeapingCSolver

class TestCSolvers(unittest.TestCase):
    """
    """

    test_model = example_models.Dimerization()
    solvers = [
        SSACSolver(model=test_model),
        ODECSolver(model=test_model),
        TauLeapingCSolver(model=test_model),
    ]
    solvers_variable = [
        SSACSolver(model=test_model, variable=True),
        ODECSolver(model=test_model, variable=True),
        TauLeapingCSolver(model=test_model, variable=True),
    ]
    
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

    @expectedFailure
    def test_solver_build(self):
        """
        Build each solver and ensure that they build properly.
        """
        expected_time = numpy.arange(100)
        print(expected_time)
        for solver in self.solvers:
            with self.subTest(solver=solver):
                results = self.test_model.run(solver=solver, number_of_trajectories=2)
                for trajectory in results:
                    self.assertTrue(trajectory["time"] == expected_time) 
