import sys
# sys.path.append('/your/gillespy2/directory/here')
if sys.version_info[0] < 3:
    raise Exception("This wrapper needs Python version 3")

import os, csv
import numpy as np
import gillespy2
from gillespy2 import TauHybridCSolver

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

# Create absolute and relative tolerance defaults
atol = 1e-12
rtol = 1e-12

# Retrieve simulation settings from file.
start, duration, steps = [0]*3
with open(SETTINGS_FILE, 'r') as settings:
    for line in settings:
        if 'start' in line: start = float(line.split(': ')[1])
        elif 'duration' in line: duration = float(line.split(': ')[1])
        elif 'steps' in line: steps = int(line.split(': ')[1])
        #elif 'absolute' in line: atol = float(line.split(': ')[1])
        #elif 'relative' in line: rtol = float(line.split(': ')[1])

#Force Continuous Species
for species in model.listOfSpecies.values():
	species.mode = 'continuous'

# Convert expected results items to species
with open(EXPECTED_RESULTS_FILE, 'r') as expected_file:
	reader = csv.reader(expected_file)
	expected_species = next(reader)[1:]
	for spec in expected_species:
		if spec in model.listOfParameters:
			p_to_s = gillespy2.Species(name=spec, initial_value=model.listOfParameters[spec].value, mode='continuous', allow_negative_populations=True)
			model.delete_parameter(spec)
			model.add_species(p_to_s)

# If no Species, plot Parameters
if not len(model.listOfSpecies):
    species_to_add = []
    for name, param in model.listOfParameters.items():
        species_to_add.append(gillespy2.Species(name=name,
                                initial_value=param.expression, mode='continuous',
                                allow_negative_populations=True))
    model.delete_all_parameters()
    model.add_species(species_to_add)
# Run simulation and store results
model.tspan = np.linspace(start, duration, steps+1)
solver = TauHybridCSolver(model=model)
results = model.run(solver=solver, show_labels=True, integrator='Radau',
integrator_options={'max_step':.25, 'atol':atol, 'rtol':rtol})

# Create headers for csv file
headers = ['time']
for species in expected_species:
    headers.append(species)

# Write results to csv file
with open(TARG_FILE, 'w+') as results_file:
    filewriter = csv.writer(results_file, delimiter=',')
    filewriter.writerow(headers)
    for timestep in range(len(results[0]['time'])):
        value_list = [results[0]['time'][timestep]]
        for species in expected_species:
            value_list.append(results[0][species][timestep])
        filewriter.writerow(value_list)

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
        print('---- Propensity_Function: ', rxn.propensity_function)

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
