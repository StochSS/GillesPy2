import unittest
import os
import tempfile
from gillespy2.core import Model, Species, Reaction, Parameter, Results, EnsembleResults

class TestResults(unittest.TestCase):

    #def name_of_test(self):
        #Put test code here 
        
    def test_to_csv_single_result_no_data(self):
        result = Results(data=None)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert not os.path.isdir(test_path)
            
    def test_to_csv_single_result_directory_exists(self):
        test_data = {'time':[0]}
        result = Results(data=test_data)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert os.path.isdir(test_path)

    def test_to_csv_single_result_file_exists(self):
        test_data = {'time':[0]}
        result = Results(data=test_data)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp+"/"+test_nametag+".csv"
            assert os.path.isfile(test_path)

    def test_to_csv_ensemble_result_no_data(self):
        result = EnsembleResults(data=None)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert not os.path.isdir(test_path)
            
    def test_to_csv_ensemble_result_directory_exists(self):
        test_data = {'time':[0]}
        result = EnsembleResults([test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp
            assert os.path.isdir(test_path)

    def test_to_csv_ensemble_result_create_one_file(self):
        test_data = {'time':[0]}
        result = EnsembleResults(data=[test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp+"/"+test_nametag+"0.csv"
            assert os.path.isfile(test_path)
        
    def test_to_csv_ensemble_result_create_multiple_files(self):
        test_data = {'time':[0]}
        result = EnsembleResults(data=[test_data,test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp = test_stamp, nametag = test_nametag, path=tempdir)
            test_path = tempdir+"/"+test_nametag+test_stamp+"/"+test_nametag+"0.csv"
            assert os.path.isfile(test_path)

if __name__ == '__main__':
    unittest.main()
