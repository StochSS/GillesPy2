#!/bin/sh
# =============================================================================
# File:    ode_wrapper_for_sts.sh
# Summary: Wrapper for running GillesPy2 ODE solvers in the SBML Test Suite
#
# This script is meant to be used by the SBML Test Suite to run GillesPy2's
# ODE solver on the SBML Test Suite test cases.  Using the GUI interface of
# the SBML Test Suite, you would open the Preferences panel, define a new
# wrapper, and in the "wrapper path" field, give the path to this file on your
# local computer.  More information about using wrappers in the SBML Test Suite
# can be found at these locations (as of 2019-09-30):
#
#  - https://github.com/sbmlteam/sbml-test-suite/wiki/Writing-wrappers-for-software
#  - https://github.com/sbmlteam/sbml-test-suite/tree/master/src/test-runner/testsuite-ui
#
# This file is a polyglot script valid in both sh and Python.  It's started
# by /bin/sh and then the quoted lines below invoke Python on this file.
# This approach came from https://unix.stackexchange.com/a/66242/141997 but
# an earlier and slightly different variant can be found in the blog posting at
# http://softwareswirl.blogspot.com/2011/06/starting-python-script-intelligently.html
# This approach is used to solve the problem that we don't know whether the
# user's "python" is Python 2 or 3, and thus we can't use an explicit path to
# python3 in the #! line, nor can we use #!/usr/bin/env python3.
# =============================================================================

''':'
# On MacOS and possibly other OSes, a non-login /bin/sh does not read the
# user's .profile, which means this script will not get the user's environment
# when it is run by SBML Test Runner and, consequently, will probably fail to
# work as expected.  The following few lines try to compensate by at least
# reading ~/.profile or ~/.bash_profile, as appropriate.
if [ -n "$BASH" ]; then
    if [ -f $HOME/.bash_profile ]; then
        . $HOME/.bash_profile
    fi
else
    if [ -f $HOME/.profile ]; then
        . $HOME/.profile
    fi
fi

# Now start Python on this file.
if type python3 >/dev/null 2>/dev/null; then
    exec python3 "$0" "$@"
else
    exec python "$0" "$@"
fi
'''

# =============================================================================
# The rest of this file is written in Python.
# =============================================================================

import sys
if sys.version_info[0] < 3:
    raise Exception("This wrapper needs Python version 3")

import os, csv
import numpy as np
import gillespy2
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver

help = \
'''This script is an application wrapper for interfacing GillesPy2's ODE solver
API to the SBML Test Suite so that the latter's Test Runner can be used to run
GillesPy2 on the semantic test cases of the SBML Test Suite.

The SBML Test Runner can supply application wrappers with the following
arguments on the command line, in any order desired:

  %d = the path to the directory containing all test cases
  %n = the current test case number (as a 5-digit number)
  %o = the directory where the CSV output file should be written
  %l = the SBML Level of the test case
  %v = the Version of SBML within the SBML Level of the test case

This wrapper script assumes they are in the same order as shown above; i.e.,
the ordering is %d %n %o %l %v.  Thus, when defining the wrapper in your copy
of the SBML Test Suite Test Runner, please define the arguments in that order.

You can optionally supply the argument -d to turn on debug printing.  This will
allow you to use the Test Runner's "View Process Output" (in the "Test" menu)
to see additional output of this wrapper when it is run by the Test Runner.
This is extremely useful for debugging and development.  The -d argument can
appear before or after the other arguments.

This script is part of GillesPy2 version {}.
For more information, please visit {}
'''.format(gillespy2.__version__, gillespy2.__url__)

if '-h' in sys.argv or '--help' in sys.argv:
    print(help)
    sys.exit()
elif len(sys.argv) < 5 or ('-d' in sys.argv and len(sys.argv) < 6):
    print('Error: missing one or more arguments.')
    print('Arguments must be TEST_PATH, CASE_NO, TARG_DIR, SBML_LVL, SBML_VER')
    print('Optional argument -d can be added to produce debugging output.')
    print('')
    print(help)
    print('Quitting.')
    sys.exit()

# Parse command line arguments and set up paths.
TEST_PATH = sys.argv[1]
CASE_NO   = sys.argv[2]
TARG_DIR  = sys.argv[3]
SBML_LVL  = sys.argv[4]
SBML_VER  = sys.argv[5]
debug     = True if '-d' in sys.argv or '--debug' in sys.argv else False

TEST_DIR = os.path.join(TEST_PATH, CASE_NO)
TEST_FILE = os.path.join(TEST_DIR, '{0}-sbml-l{1}v{2}.xml'.format(CASE_NO, SBML_LVL, SBML_VER))
SETTINGS_FILE = os.path.join(TEST_DIR, '{0}-settings.txt'.format(CASE_NO))
TARG_FILE = os.path.join(TARG_DIR, '{0}-results.csv'.format(CASE_NO))
EXPECTED_RESULTS_FILE = os.path.join(TEST_DIR, '{0}-results.csv'.format(CASE_NO))

# Read in the SBML model.
model, errors = gillespy2.import_SBML(TEST_FILE)

# Retrieve simulation times from settings file.
start, duration, steps = [0]*3
with open(SETTINGS_FILE, 'r') as settings:
    for line in settings:
        if 'start' in line: start = float(line.split(': ')[1])
        elif 'duration' in line: duration = float(line.split(': ')[1])
        elif 'steps' in line: steps = int(line.split(': ')[1])

# Run simulation and store results
solver = BasicODESolver()
model.tspan = np.linspace(start, duration, steps+1)
results = model.run(solver=solver, show_labels=False)

# Create headers for csv file
headers = ['time']
for species in model.listOfSpecies:
    headers.append(species)

# Write results to csv file
with open(TARG_FILE, 'w+') as results_file:
    filewriter = csv.writer(results_file, delimiter=',')
    filewriter.writerow(headers)
    for row in results[0]:
        filewriter.writerow(row)

if debug:
    print('SBML FILE:')
    with open(TEST_FILE, 'r') as sbml_file:
        print(sbml_file.read())

    print('Solver: ', solver.name)
    print('Model: ', model.name)

    print('\nStart: ', start)
    print('Duration: ', duration)
    print('Steps: ', steps)

    print('\nParameters:')
    for name, param in model.listOfParameters.items():
        print('-- {0}: {1}'.format(name, param.expression))
    print('\nSpecies:')
    for name, spec in model.listOfSpecies.items():
        print('-- {0}: {1}'.format(name, spec.initial_value))
    print('\nReactions:')
    for name, rxn in model.listOfReactions.items():
        print('-- ', name)
        print('---- Reactants:')
        for rct, cnt in rxn.reactants.items():
            print('------ {0}: {1}'.format(str(rct), cnt))
        print('---- Products:')
        for prd, cnt in rxn.products.items():
            print('------ {0}: {1}'.format(str(prd), cnt))

    print('Simulation Results:')
    with open(TARG_FILE, 'r') as results_file:
        for row in results_file:
            print(row)
    print('Expected Results:')
    with open(EXPECTED_RESULTS_FILE, 'r') as expected_file:
        for row in expected_file:
            print(row)



# Please leave the following lines for Emacs users
# ......................................................................
# Local Variables:
# mode: python
# python-indent-offset: 4
# End:
