# encoding: utf-8
from gillespy2.core import GillesPySolver, gillespyError, log
import numpy as np
import random
cimport numpy as np
cimport cython
from cpython.exc cimport PyErr_CheckSignals
from libc.stdlib cimport malloc, free
cimport libc.math as math
import re

cdef struct Equation:
    int length
    int *terms
    int parameters_length
    double *parameters
    

cdef struct CythonReaction:
    int *affected_reactions
    Equation propensity_function

cdef enum Operators:
    add = 1
    sub = 2
    mul = 3
    div = 4
    exp = 5
    parameter = 6
    
DEF MAX_STACK_SIZE=100
cdef double operand_stack[MAX_STACK_SIZE]
cdef int operator_stack[MAX_STACK_SIZE]


cdef void simulate_trajectory(np.ndarray[np.float64_t, ndim=2] trajectory, CythonReaction *reactions, int number_reactions, np.ndarray[np.float64_t, ndim=2] species_changes, int seed):
    cdef int i,j
    cdef double current_time = 0
    cdef int number_entries = 0
    cdef int rc = 0
    cdef np.ndarray[np.float64_t, ndim=1] current_state = np.zeros((trajectory.shape[1]-1))
    cdef double *propensities = <double*> malloc(number_reactions * sizeof(double))
    np.copyto(current_state, trajectory[0,1:])
    for i in range(number_reactions):
        propensities[i] = evaluate_prefix(reactions[i].propensity_function, (current_state))
    cdef double propensity_sum, cumulative_sum
    while number_entries < trajectory.shape[0]:
        PyErr_CheckSignals()
        propensity_sum = 0
        for i in range(number_reactions):
            propensity_sum += propensities[i]
        if propensity_sum <= 0:
            trajectory[number_entries:,1:] = current_state
            break
        if seed >= 0:
            random.seed(seed)
        cumulative_sum = random.random() * propensity_sum
        current_time -= math.log(random.random()) / propensity_sum
        while number_entries < trajectory.shape[0] and trajectory[number_entries, 0] <= current_time:
            trajectory[number_entries, 1:] = current_state
            number_entries += 1
        for i in range(number_reactions):
            cumulative_sum -= propensities[i]
            if cumulative_sum <= 0:
                current_state += species_changes[i]
                for j in range(number_reactions):
                    propensities[j] = evaluate_prefix(reactions[j].propensity_function, (current_state))
                break
    free(propensities)

#Evaluates an equation in Polish notation from left to right
@cython.boundscheck(False)
cdef double evaluate_prefix(Equation eqn, np.ndarray[np.float64_t, ndim=1] state):
    cdef int i = 0
    cdef int operands = 0
    cdef int operators = 0
    cdef double current_operand = 0
    cdef int pending_operand = 0
    for i in range(eqn.length):
        if eqn.terms[i] >= 0 and eqn.terms[i] < Operators.parameter:
            operator_stack[operators] = eqn.terms[i]
            operators += 1
            pending_operand = 0
        else:
            if eqn.terms[i] >= Operators.parameter:
                current_operand = eqn.parameters[eqn.terms[i]-Operators.parameter]
            elif eqn.terms[i] < 0:
                current_operand = state[-(eqn.terms[i] + 1)]
            if pending_operand:
                while operands > 0:
                    operands -= 1
                    operators -= 1
                    if operator_stack[operators] == Operators.add:
                        current_operand = operand_stack[operands] + current_operand
                    elif operator_stack[operators] == Operators.sub:
                        current_operand = operand_stack[operands] - current_operand
                    elif operator_stack[operators] == Operators.mul:
                        current_operand = operand_stack[operands] * current_operand
                    elif operator_stack[operators] == Operators.div:
                        current_operand = operand_stack[operands] / current_operand
                    elif operator_stack[operators] == Operators.exp:
                        current_operand = math.pow(operand_stack[operands], current_operand)
            pending_operand = 1
            operand_stack[operands] = current_operand
            operands += 1
    return operand_stack[0]

operator_precedence = {
    '^' : (4, 'R'),
    '*' : (3, 'L'),
    '/' : (3, 'L'),
    '+' : (2, 'L'),
    '-' : (2, 'L')
}
def convert_infix_prefix(equation):
    #stack for final result
    output = []
    #stack for operators parsed
    operators = []
    while len(equation) > 0:
        token = equation.pop()
        if token in operator_precedence:
            while len(operators) > 0 and operators[0] != ')' and operator_precedence[operators[0]][1] == 'L' and operator_precedence[operators[0]][0] > operator_precedence[token][0]:
                output.insert(0, operators.pop(0))
            operators.insert(0,token)
        elif token == ')':
            operators.insert(0,token)
        elif token == '(':
            while len(operators) > 0 and operators[0] != ')':
                output.insert(0, operators.pop(0))
            if len(operators) > 0 and operators[0] == ')':
                operators.pop(0)
        else:
            output.insert(0, token)
    while len(operators) > 0:
        output.insert(0, operators.pop(0))
    return output
                
    
class CythonSSASolver(GillesPySolver):
    '''
    This solver was deprecated in release version 1.3 of GillesPy2.
    '''
    name = "CythonSSASolver"
    interrupted = False
    rc = 0

    def __init__(self):
        name = "CythonSSASolver"
        interrupted = False
        rc = 0

    #@cython.boundscheck(False)
    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, profile=False, show_labels=True, **kwargs):

        import signal
        def timed_out(signum, frame):
            print('signal raised')
            self.rc = 33
            self.interrupted = True
            raise gillespyError.SimulationTimeoutError()

        signal.signal(signal.SIGALRM, timed_out)

        if not isinstance(self, CythonSSASolver):
            self = CythonSSASolver()

        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))

        self.simulation_data = []
        #convert dictionary of species to species array
        species = list(model.listOfSpecies.keys())
        cdef int number_species = len(species)

        if seed is not None:
            if not isinstance(seed, int):
                seed = int(seed)
            if seed < 0:
                raise gillespyError.ModelError('seed must be a positive integer')
        else:
            seed = -1
        cdef int seed_arg = seed
        #set timespan for simulation(s)
        timeline = np.linspace(0,t, int(round(t /increment+1)))
        #allocate memory for trajectories
        cdef np.ndarray[np.float64_t, ndim=3] trajectories = np.zeros((number_of_trajectories, timeline.size, number_species + 1))
        trajectories[:,:,0] = timeline
        cdef int i = 0, j
        for i in range(number_species):
            trajectories[:,:,i+1] = model.listOfSpecies[species[i]].initial_value
    #convert dictionary of reactions to reactions array
        cdef int number_reactions = len(list(model.listOfReactions.keys()))
        cdef CythonReaction *reactions = <CythonReaction*> malloc(number_reactions * sizeof(CythonReaction))
        i = 0
        reaction_names = list(model.listOfReactions.keys())
        for i in range(number_reactions):
            reactions[i].affected_reactions = <int*> malloc(number_reactions*sizeof(int))

    #convert propensity functions now
        #create dictionary of all constant parameters for propensity evaluation
        parameters = {'vol' : model.volume}
        for paramName, param in model.listOfParameters.items():
            parameters[paramName] = param.value
        propensity_functions = [" "+r.propensity_function.replace(' ','') for r in model.listOfReactions.values()]
        #get all numeric constants from propensity functions
	# TODO THIS BLOCK OF CODE CONTAINS AN RE MATCHING ERROR FOR PROP FUNCTIONS CONTAINING SUBTRACTION
        numbers = re.compile('[^\[\]\w](\-?\d+\.?\d*)')
        constants = []
        for fun in propensity_functions:
            matches = numbers.findall(fun)
            for match in matches:
                if match not in constants:
                    constants.append(match[:])
        for c in constants:
            parameters[c] = float(c)
        paramNames = list(parameters.keys())
	# TODO END ERROR BLOCK
        #sort parameter names by longest first to prevent issues of one param containing a shorter param in name
        paramNames.sort(key = lambda x: -len(x))
        cdef double *cParameters = <double*> malloc(len(paramNames)*sizeof(double))
        for i in range(len(paramNames)):
            cParameters[i] = parameters[paramNames[i]]
        #create regex for parsing propensity function tokens
        equation_terms = re.compile('[\+\-\*/\^\(\)]|x\d+|y\d+')
        cdef np.ndarray[np.float64_t, ndim=2] species_changes = np.zeros((number_reactions, number_species))
        #pre-evaluate propensity equations from strings:
        for i in range(number_reactions):
            reactions[i].propensity_function.parameters = cParameters
            reactions[i].propensity_function.parameters_length = len(paramNames)
            reaction = model.listOfReactions[reaction_names[i]]
            #replace all references to parameters/constants
            for j in range(len(paramNames)):
                propensity_functions[i] = propensity_functions[i].replace(paramNames[j],'y{0}'.format(j))
            #replace all references to species with array indices
            spec_index = {}
            spec_list = []
            for j in range(number_species):
                spec = model.listOfSpecies[species[j]]
                spec_index[spec] = j
                spec_list.append(spec)
            spec_list.sort(key=lambda x: -len(str(x)))
            for spec in spec_list:
                species_changes[i][spec_index[spec]] = reaction.products.get(spec,0) - reaction.reactants.get(spec, 0)
                propensity_functions[i] = propensity_functions[i].replace(str(spec), 'x{0}'.format(spec_index[spec]))
            prefix_eqn = convert_infix_prefix(equation_terms.findall(propensity_functions[i]))
            reactions[i].propensity_function.length = len(prefix_eqn)
            reactions[i].propensity_function.terms = <int*> malloc(reactions[i].propensity_function.length * sizeof(int))
            for j in range(reactions[i].propensity_function.length):
                if prefix_eqn[j] == '+':
                    reactions[i].propensity_function.terms[j] = Operators.add
                elif prefix_eqn[j] == '-':
                    reactions[i].propensity_function.terms[j] = Operators.sub
                elif prefix_eqn[j] == '*':
                    reactions[i].propensity_function.terms[j] = Operators.mul
                elif prefix_eqn[j] == '/':
                    reactions[i].propensity_function.terms[j] = Operators.div
                elif prefix_eqn[j] == '^':
                    reactions[i].propensity_function.terms[j] = Operators.exp
                elif prefix_eqn[j][0] == 'x':
                    reactions[i].propensity_function.terms[j] = -1 - int(prefix_eqn[j][1:])
                elif prefix_eqn[j][0] == 'y':
                    reactions[i].propensity_function.terms[j] = Operators.parameter + int(prefix_eqn[j][1:])
                else:
                    print("Error: Unrecognized term: {0}".format(prefix_eqn[j]))
        #begin simulating each trajectory
        cdef double current_time, propensity_sum, cumulative_sum
        cdef np.ndarray[np.float64_t, ndim=1] current_state = np.zeros((number_species))
        cdef int number_threads = 4
        for i in range(number_of_trajectories):
            simulate_trajectory(trajectories[i], reactions, number_reactions, species_changes, seed_arg)
            print(self.rc)
            #assemble complete simulation data in format specified
            if show_labels:
                data = {'time' : timeline}
                for j in range(number_species):
                    data[species[j]] = trajectories[i,:,j+1]
                self.simulation_data.append(data)
            else:
                self.simulation_data = trajectories
        #clean up
        for i in range(number_reactions):
            free(reactions[i].affected_reactions)
            free(reactions[i].propensity_function.terms)
        free(reactions)
        free(cParameters)
        return self.simulation_data, self.rc
        
