import os, sys, csv
import matplotlib.pyplot as plt
import numpy as np
np.set_printoptions(suppress=True)


'''
This script will scan the local directory for sbml test cases (stored in folders named after the test case, and will also scan a results directory.  It compares the expected test case results to the generated results in the results folder.
'''

print('Test Numbers for completed tests:\n')
for name in os.listdir('./results'):
    if 'results' in name:
        print(name.split('-')[0])
#Arguments to call should be: TEST_NO

#TEST_NO = sys.argv[1]
TEST_NO = input('\nCompare which results? Enter Test_number: ')
print('test')
SOLVER_NO = input('\nCompare which results? Enter Solver_number: \
                   \n 1. Auto \
                   \n 2. Hybrid \n\n')
if SOLVER_NO.lower() is '1' or 'auto':
    TARG_DIR = os.path.join(os.getcwd(), 'auto_results')
elif SOLVER_NO.lower() is '2' or 'hybrid':
    TARG_DIR = os.path.join(os.getcwd(), 'hybrid_results')
    
EXPECTED_DIR = os.path.join(os.getcwd(), TEST_NO)
EXPECTED_RESULTS = os.path.join(EXPECTED_DIR, '{0}-results.csv'.format(TEST_NO))
#TARG_DIR = os.path.join(os.getcwd(), 'results')
SIM_RESULTS = os.path.join(TARG_DIR, '{0}-results.csv'.format(TEST_NO))
TARG_FILE = os.path.join(TARG_DIR, '{0}-scores.csv'.format(TEST_NO))

expected_data = []

with open(EXPECTED_RESULTS, 'r') as expected:
    reader = csv.reader(expected, delimiter=',')
    dtype= np.float64
    headers = next(reader)
    #t = []
    #spec_names = []
    #for item in headers[1:]:
        #spec_names.append(item)
    #print(spec_names)
    for i, row in enumerate(reader):
        if any(x.strip() for x in row):
            #t.append(row[0])
            #for item in row[1:]:
                #plt.plot(row[0], item)
            expected_data.append(np.array(row, dtype=dtype))

expected_data = np.array(expected_data)


sim_data = []
with open(SIM_RESULTS, 'r') as sim:
    reader = csv.reader(sim, delimiter=',')
    dtype= np.float64
    next(reader)
    for row in reader:
        sim_data.append(np.array(row, dtype=dtype))

sim_data = np.array(sim_data)
print('\nExpected:')
print(headers)
print(expected_data)

print('\nSimulation:')
print(headers)
print(sim_data)

results = expected_data
result = 'PASS'
for i in range(len(expected_data)):
    for j in range(len(expected_data[i])):
        if j == 0:
            results[i][j] = sim_data[i][j] # copy time
        elif expected_data[i][j] == 0 and sim_data[i][j] == 0:
            results[i][j] = 0
        else:
            results[i][j] = abs((expected_data[i][j] / sim_data[i][j])-1) * 100
            if results[i][j] > 3:
                result = 'FAIL'

print('\nPercent difference')
print(headers)
print(results)
print()
print(result)

with open(TARG_FILE, 'w+') as scores:
    filewriter = csv.writer(scores, delimiter=',')
    filewriter.writerow(headers)
    for row in results:
        filewriter.writerow(row)
