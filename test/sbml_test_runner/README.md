# Running the SBML Compliance Tests  
  
In order to verify SBML compliance with the GillesPy2 package, follow these steps.  
  
## 1. Install the [SBML Test Suite](http://sbml.org/Software/SBML_Test_Suite)  
  
## 2. Configure wrappers. 
    While running the SBML Test Suite, go to File->Options/Wrappers  
    The wrappers for gillespy2 solvers are contained [here](/test/sbml_test_runner/wrappers): /test/sbml_test_runner/wrappers  
    Provide the path for the desired wrapper along with a name.
    Be sure to configure all arguments to the wrapper. The 'arguments to wrapper' field should contain '%d %n %o %l %v'
    If you would like additional information to be printed to stdout while running the tests, also include '-d' at the end of 'arguments to wrapper'
    You can provide any output directory you desire.

## 3. In the file for the wrapper you are using, be sure to append your GillesPy2 directory to the system path

## 4. Save, and begin running test cases.  Additional documentation can be found at [the SBML site](http://sbml.org/Software/SBML_Test_Suite) or [the SBML github repo](https://github.com/sbmlteam/sbml-test-suite/)
