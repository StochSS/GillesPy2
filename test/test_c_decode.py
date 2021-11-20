from datetime import time
import unittest
import numpy as np
import io
from gillespy2.solvers import ODECSolver, SSACSolver, TauLeapingCSolver, TauHybridCSolver
from example_models import Degradation
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder, IterativeSimDecoder

class TestCDecode(unittest.TestCase):
    model = Degradation()
    solvers = [
        # Solver     List of trajectory numbers to try
        (ODECSolver,        [1]),
        (SSACSolver,        [1, 3]),
        (TauLeapingCSolver, [1, 3]),
        (TauHybridCSolver,  [1, 3]),
    ]
    decoders = [BasicSimDecoder, IterativeSimDecoder]

    def setUp(self):
        self.model.timespan(np.linspace(0, 300, 301))
        self.solvers = [(solver(model=self.model), trajectory_numbers) for solver, trajectory_numbers in self.solvers]

    def test_run_output(self):
        expected_tspan = self.model.tspan
        expected_init = self.model.get_species("A").initial_value
        for solver, trajectory_counts in self.solvers:
            for number_of_trajectories in trajectory_counts:
                with self.subTest("Processing simulation output for each solver with different trajectory sizes",
                                  number_of_trajectories=number_of_trajectories,
                                  solver=solver):
                    results = self.model.run(solver=solver, number_of_trajectories=number_of_trajectories, seed=1024)
                    for result in results:
                        tspan, first_value, last_value = result["time"], result["A"][0], result["A"][-1]
                        result_diff = np.concatenate([np.array([first_value]), result["A"][1:]])
                        result_diff = result["A"] - result_diff
                        self.assertTrue(np.allclose(tspan, expected_tspan), msg="C++ Simulation output contains unexpected timeline values"
                                                                        f"\n  Received timeline: {tspan}\n  Expected timeline: {expected_tspan}")
                        self.assertEqual(first_value, expected_init, msg=f"C++ Simulation output begins with unexpected value: {first_value}")
                        self.assertAlmostEqual(last_value, 0, places=3, msg=f"C++ Simulation output converges on an unexpectedly large value: {last_value}")

    def test_c_decoder(self):
        """
        Ensures that, given a certain input from stdout, the results will be
            properly formatted by the sim decoders.
        """
        example_input = """
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,
        0,1.0,2.0,3.0,1,2.0,4.0,6.0,2,3.0,6.0,9.0,3,4.0,8.0,12.0,3
        """.strip()
        expected_result = np.array([
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ],
            [ [0, 1.0,2.0,3.0], [1, 2.0,4.0,6.0], [2, 3.0,6.0,9.0], [3, 4.0,8.0,12.0] ]
        ])

        for decoder in TestCDecode.decoders:
            mock_stdout = io.BytesIO(bytes(example_input, "ascii"))
            mock_stdout = io.BufferedReader(mock_stdout)
            with self.subTest("Processing mock data output with C decoders", decoder=decoder):
                reader = decoder(np.zeros((2, 4, 4)))
                reader.read(mock_stdout)
                result, time_stopped = reader.get_output()

                self.assertEqual(time_stopped, 3, f"Stop time returned by decoder is incorrect (expected {3}, got {time_stopped})")
                # Test both the passed-in trajectory list and the returned trajectories.
                self.assertTrue(np.all(result == expected_result), "Trajectory data returned by decoder is incorrect")
