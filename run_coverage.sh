#!/bin/bash
coverage run test/run_tests.py -m release
coverage html
