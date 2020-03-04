# CONTRIBUTING


If you are interested in contributing to GillesPy2, please follow these guidelines:  

### 1. Making Changes:  
a) Create a fork from this repository ('https://github.com/GillesPy2/GillesPy2')  
b) Create a new branch.  
c) Please use a meaningful name, and include the issue # if applicable.  
d) For bug fixes, branch from 'master'  
e) For new features, branch from 'develop'  
e) Be sure to document your code  
f) Create unit tests to cover changes.  Unit tests should be placed in the [test](/test) subdirectory.
g) Update the UML class model and UML class diagram screenshot located in the [docs/getting_started/uml_class_model](/docs/getting_started/uml_class_model) subdirectory. The model was built with Pynsource, an open-source UML tool for Python.
  
### 2. Submitting a Pull Request:  
a) If changes are bug/hotfix, make a pull request to 'master'.  
b) If other changes, or new features are being added, make a pull request to 'develop'.  
c) Include a list of changes in the PR description.  
d) Provide a usage guide/how-to for new features in the PR description.  
  
***Do NOT merge your own pull request.  Once the request is submitted, the code must be reviewed and merged by someone else.***  
  
### 3. Merging a Pull Request:
a) Verify correct merge destination ('master' for hotfix/bugs, 'develop' for features/changes).  
b) Review code for logic, consistency, and documentation.  
c) Verify that unit tests are provided for the new code, and that they accurately test the new feature/changes.  
d) Check coverage by running ***./run_coverage.sh*** from the project root directory. The coverage results are stored in the htmlcov subdirectory, which can be compared with the files changed (viewable through github) to ensure that the changes made have been covered by proper testing.  
e) Make sure that all unit tests in Travis CI pass before merging the changes.  
f) Merge!  

