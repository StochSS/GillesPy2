#! /usr/bin/python3.6

import os, sys, csv
import numpy as np
sys.path.append('/home/smatthe2/GillesPy2')
import gillespy2
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver

# Parse Args: TEST_PATH, CASE_NO, TARG_DIR, SBML_LVL, SBML_VER

assert len(sys.argv) > 5, 'Did not provide all arguments [5]: TEST_PATH, CASE_NO, TARG_DIR, SBML_LVL, SBML_VER, [-d]'

#Parse command line arguments and set up paths
TEST_PATH = sys.argv[1]
CASE_NO = sys.argv[2]
TARG_DIR = sys.argv[3]
SBML_LVL = sys.argv[4]
SBML_VER = sys.argv[5]
if len(sys.argv) > 6 and sys.argv[6] == '-d': debug = True
else: debug = False
TEST_DIR = os.path.join(TEST_PATH, CASE_NO)
TEST_FILE = os.path.join(TEST_DIR, '{0}-sbml-l{1}v{2}.xml'.format(CASE_NO, SBML_LVL, SBML_VER))
SETTINGS_FILE = os.path.join(TEST_DIR, '{0}-settings.txt'.format(CASE_NO))
TARG_FILE = os.path.join(TARG_DIR, '{0}-results.csv'.format(CASE_NO))
EXPECTED_RESULTS_FILE = os.path.join(TEST_DIR, '{0}-results.csv'.format(CASE_NO))

solver = BasicODESolver()
model, errors = gillespy2.import_SBML(TEST_FILE)

# retrieve simulation times from settings file
start, duration, steps = [0]*3
with open(SETTINGS_FILE, 'r') as settings:
    for line in settings:
        if 'start' in line: start = float(line.split(': ')[1])
        elif 'duration' in line: duration = float(line.split(': ')[1])
        elif 'steps' in line: steps = int(line.split(': ')[1])

# Run simulation and store results
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


