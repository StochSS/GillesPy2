#!/bin/bash
#pip3 install  coverage
#pip3 install  coverage-badge

coverage run --source=gillespy2 --omit=gillespy2/solvers/stochkit/* test/run_tests.py -m develop
coverage html
coverage-badge -fo coverage.svg
