import unittest
from gillespy2.core import Model, Species, Reaction, Parameter, Results, EnsembleResults
from gillespy2.core.gillespyError import *
import numpy as np
import os
from datetime import datetime

class TestResults(unittest.TestCase):

    #def name_of_test(self):
        #Put test code here 
        
    def test_to_csv_single(self):
        result = new Result()
        now = datetime.now()
        timestamp=datetime.timestamp(now)
        test_nametag = "test"
        with TemporaryDirectory() as tempdir:
            result.to_csv(nametag = test_nametag, path=tempdir)
            test_path = os.path.join(tempdir,test_nametag+str(timestamp))
            assert os.path.isdir(test_path)
            assert os.path.isfile(test_path + "/" + test_nametag)
            assert test_nametag == "fish"
            
    def test_to_csv_ensemble(self):
        result = new EnsembleResults()
        
        

if __name__ == '__main__':
    unittest.main()
