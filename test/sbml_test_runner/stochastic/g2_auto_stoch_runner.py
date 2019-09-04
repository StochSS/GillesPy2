#! /usr/bin/python3.6

import os, sys, csv
import numpy as np
sys.path.append('/home/smatthe2/GillesPy2')
import gillespy2
from gillespy2.solvers.numpy.ssa_solver import NumPySSASolver
'''
This script will scan all subdirectories it shares a directory with for sbml case tests.  The user will be prompted to select from all discovereted test cases, and then to select a number of trajectories to run.  The results will be stored in the results subdirectory in csv format for future comparison.
'''

def run_simulation(test_file, targ_file, expected_results_file, num_sims, timeout ):
    print('now simulating ', test_file)
    model, errors = gillespy2.import_SBML(test_file)
    solver = NumPySSASolver()

    # retrieve simulation times from settings file
    start, duration, steps = [0]*3
    with open(settings_file, 'r') as settings:
        for line in settings:
            if 'start' in line: start = float(line.split(': ')[1])
            elif 'duration' in line: duration = float(line.split(': ')[1])
            elif 'steps' in line: steps = int(line.split(': ')[1])

    # Run simulation and store results
    model.tspan = np.linspace(start, duration, steps+1)

    if debug:
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

    results = model.run(solver=solver, time_out=TIME_OUT, show_labels=False, number_of_trajectories=NUMBER_OF_SIMS)
    averages = np.mean(results, axis=0)
    std = np.std(results, axis=0)
    avg_only = averages[:, 1:]
    std_only = std[:, 1:]
    combined = np.concatenate((avg_only, std_only), axis=1)
    times = np.array(results[0, :, 0])
    data = np.insert(combined, (0), times, axis=1)
    print('data: ', data.shape)
    print(data)

    # Create headers for csv file
    headers = ['time']
    for species in model.listOfSpecies:
        headers.append('{0}-mean'.format(species))
    for species in model.listOfSpecies:
        headers.append('{0}-sd'.format(species))

    # Write results to csv file
    with open(targ_file, 'w+') as results_file:
        filewriter = csv.writer(results_file, delimiter=',')
        filewriter.writerow(headers)
        for row in data:
            filewriter.writerow(row)

    if debug:
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
        with open(targ_file, 'r') as results_file:
            for row in results_file:
                print(row)
        print('Expected Results:')
        with open(expected_results_file, 'r') as expected_file:
            for row in expected_file:
                print(row)



print('LIST OF TEST CASES\n')
for name in os.listdir():
    if os.path.isdir(name) and '0' in name:
        print(name)
TEST_PATH = '.'
CASE_NO = input('\nCase Number, or type "ALL" to run all tests: ')
TIME_OUT = int(input('\nPlease input a simulation timeout duration (in seconds) per trajectory.  Use "0" to disable timeout: '))
TARG_DIR = './results/'
SBML_LVL = '3'
SBML_VER = '1'
NUMBER_OF_SIMS = int(input('Number of Simulations to run: '))

debug = True

if CASE_NO.lower() == 'all':
    print('all selected')
    for case_no in os.listdir():
        if os.path.isdir(case_no) and '0' in case_no:
            print('processing ', case_no)
            test_dir = os.path.join(TEST_PATH, case_no)
            test_file = os.path.join(test_dir, '{0}-sbml-l{1}v{2}.xml'.format(case_no, SBML_LVL, SBML_VER))
            settings_file = os.path.join(test_dir, '{0}-settings.txt'.format(case_no))
            targ_file = os.path.join(TARG_DIR, '{0}-results.csv'.format(case_no))
            expected_results_file = os.path.join(test_dir, '{0}-results.csv'.format(case_no))
            try:
                run_simulation(test_file, targ_file, expected_results_file, NUMBER_OF_SIMS, TIME_OUT)
            except:
                print('could not complete test case for ', test_file)
else:
    TEST_DIR = os.path.join(TEST_PATH, CASE_NO)
    TEST_FILE = os.path.join(TEST_DIR, '{0}-sbml-l{1}v{2}.xml'.format(CASE_NO, SBML_LVL, SBML_VER))
    SETTINGS_FILE = os.path.join(TEST_DIR, '{0}-settings.txt'.format(CASE_NO))
    TARG_FILE = os.path.join(TARG_DIR, '{0}-results.csv'.format(CASE_NO))
    EXPECTED_RESULTS_FILE = os.path.join(TEST_DIR, '{0}-results.csv'.format(CASE_NO))
    run_simulation(TEST_FILE, TARG_FILE, EXPECTED_RESULTS_FILE, NUMBER_OF_SIMS, TIME_OUT)
    

