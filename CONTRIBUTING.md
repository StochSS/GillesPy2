# CONTRIBUTING

By contributing to GillesPy2, you agree to the [Code of Conduct](/CODE_OF_CONDUCT.md) and that the contents of this project are distributed under the [GNU GPL v3](/LICENSE). Please review both of these documents. If after reading you are interested in contributing, please follow these guidelines:  

### 1. Making Changes:  
a) Create a fork from this repository ('https://github.com/GillesPy2/GillesPy2')  
b) Create a new branch.  
c) Please use a meaningful name, and include the issue # if applicable.  
d) For bug fixes, branch from 'master'  
e) For new features, branch from 'develop'  
f) Be sure to document your code  
g) If you are creating a new file, prepend the [COPYRIGHT](/COPYRIGHT) notice to the beginning of the file.  
h) Create unit tests to cover changes. Unit tests should be placed in the [test](/test) subdirectory.  
i) Update [UML_CLASS_DIAGRAM.md](UML_CLASS_DIAGRAM.md) in the root directory, the [UML class diagram](.graphics/gillespy2-UML-class-diagram.png) in the [.graphics](.graphics) subdirectory, and the [UML class model](docs/getting_started/basic_usage/gillespy2-UML-class-model.pyns) in the [docs/getting_started/basic_usage/](docs/getting_started/basic_usage/) subdirectory to reflect your changes. These files were created using [Pynsource](https://pynsource.com/), an open-source UML tool for Python.
  
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
