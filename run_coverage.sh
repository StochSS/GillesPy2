#!/bin/bash
#pip3 install  coverage
#pip3 install python-libsbml

coverage run --source=gillespy2 --omit=*gillespy2/solvers/stochkit/*,*gillespy2/solvers/cython/* test/run_tests.py -m develop
coverage html
