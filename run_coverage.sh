#!/bin/bash
coverage run test/run_tests.py -m develop
coverage html
coverage-badge -fo coverage.svg
