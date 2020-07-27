# Release Checklist


### Github/Code

1. Ensure requirements.txt (dependencies) and our acknowledgments are up to date

2. Be sure to mark code as deprecated/remove as necessary 

3. Verify that user documentation matches the current release

4. Verify all testing procedures 

5. Update authors/maintainers list and acknowledgments (in README.md as well as version.py!)

6. Update all instances of the version number (including version.py)

7. Push changes to 'main'




#############################################
## This is done automatically though Travis CI, ensure release is tagged on github
#############################################
### PyPI

1. After following the above, ensure that setup.py will be updated correctly per changes in version.py, 
   ensure any other necessary changes to setup.py
   
2. Build wheel and source file for upload to PyPI

3. Follow instructions on link below to TEST our PyPI package (refer to #6 on this list after doing so)

4. Upload to PyPI

5. Ensure pip install works properly 

6. Check PyPI listing page to match that README and other documentation display properly


[For documentation on testing PyPI, click here](https://packaging.python.org/guides/using-testpypi/#using-test-pypi)

[For complete documentation on uploading to PyPI, including code, refer to this link](https://packaging.python.org/guides/distributing-packages-using-setuptools/)
#############################################



### conda-forge

1. Ensure you are entered as a maintainer in the meta.yaml file of our conda-forge fork

2. Create fork of our feedstock at the conda-forge feedstocks github

3. In fork, update all necessary elements of meta.yaml. This includes:
   authors, version numbers, URL and sha256, dependencies
   
4. Bump build number AS NEEDED (click on first link for appropriate reasons to do so)
   
4. Commit changes and submit pull request to feedstock location

[For complete documentation on updating conda-forge feedstock, including code, refer to this link](https://conda-forge.org/docs/maintainer/updating_pkgs.html#)

[For complete information on various attributes of the meta.yaml file, refer to this link](https://conda-forge.org/docs/maintainer/adding_pkgs.html)
