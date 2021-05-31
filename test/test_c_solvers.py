import unittest
import numpy
import io
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder

class TestCSolvers(unittest.TestCase):
    """
    """
    
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

