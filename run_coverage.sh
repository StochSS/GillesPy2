#!/bin/bash
pip3 install --user coverage-badge

coverage run test/run_tests.py -m develop
coverage html
coverage-badge -fo coverage.svg
