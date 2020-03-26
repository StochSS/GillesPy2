#!/bin/bash
#pip3 install  coverage
#pip3 install  coverage-badge
#pip3 install python-libsbml

coverage run --source=gillespy2 --omit=gillespy2/solvers/stochkit/* test/run_tests.py -m develop
coverage html
