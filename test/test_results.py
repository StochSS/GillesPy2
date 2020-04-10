import unittest
import os
import tempfile
from gillespy2.core import Model, Species, Reaction, Parameter, Results, EnsembleResults

class TestResults(unittest.TestCase):
   
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

    def test_to_csv_single_result_no_stamp(self):
        test_data = {'time':[0]}
        result = Results(data=test_data)
        test_nametag = "test_nametag"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(nametag=test_nametag, path=tempdir)
            assert len(os.listdir(tempdir)) is not 0

    def test_to_csv_single_result_no_nametag(self):
        test_data = {'time':[0]}
        test_model = Model('test_model')
        result = Results(data=test_data,model=test_model)
        result.solver_name = 'test_solver'
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(stamp=test_stamp, path=tempdir)
            assert len(os.listdir(tempdir)) is not 0

    def test_to_csv_single_result_no_path(self):
        test_data = {'time':[0]}
        result = Results(data=test_data)
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            result.to_csv(stamp=test_stamp, nametag=test_nametag)
            assert len(os.listdir(tempdir)) is not 0

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

    def test_to_csv_ensemble_result_no_stamp(self):
        test_data = {'time':[0]}
        result = EnsembleResults(data=[test_data,test_data])
        test_nametag = "test_nametag"

        with tempfile.TemporaryDirectory() as tempdir:
            result.to_csv(nametag=test_nametag, path=tempdir)
            assert len(os.listdir(tempdir)) is not 0

    def test_to_csv_ensemble_result_no_nametag(self):
        test_data = {'time':[0]}
        test_model = Model('test_model')
        result = Results(data=test_data,model=test_model)
        ensemble_result = EnsembleResults(data=[result,result])
        result.solver_name = 'test_solver'
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            ensemble_result.to_csv(stamp=test_stamp, path=tempdir,nametag=None)
            assert len(os.listdir(tempdir)) is not 0

    def test_to_csv_ensemble_result_no_path(self):
        test_data = {'time':[0]}
        result = EnsembleResults(data=[test_data,test_data])
        test_nametag = "test_nametag"
        test_stamp = "test_stamp"

        with tempfile.TemporaryDirectory() as tempdir:
            os.chdir(tempdir)
            result.to_csv(stamp=test_stamp, nametag=test_nametag)
            assert len(os.listdir(tempdir)) is not 0

if __name__ == '__main__':
    unittest.main()
