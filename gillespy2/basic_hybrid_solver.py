from .gillespySolver import GillesPySolver
import random
import math
import numpy
from scipy.integrate import ode
import scipy.integrate

eval_globals ={}

class BasicHybridSolver(GillesPySolver):
    name = "Basic Hybrid Solver"

    def __init__(self, debug=False):
        self.debug = debug
    
    
    @staticmethod
    # Create the matrix equation for integration
    def f(t, y, curr_state, species, parameters, reactions, rate_rules, eval_globals):
        #curr_state = {'vol': 1}
        curr_state['t'] = t
        state_change = []
        #for i, s in enumerate(species):
        #    curr_state[s] = pops[i]
        #for p in parameters:
        #    curr_state[p] = parameters[p].value
        for i, r in enumerate(reactions):
            state_change.append(eval(reactions[r].propensity_function, eval_globals, curr_state))
        for i, rr in enumerate(rate_rules):
            state_change.append(eval(rate_rules[rr].expression,  eval_globals, curr_state))

        return state_change

#    @staticmethod
#    def get_reaction(curr_state, y0, model, curr_time, save_time, eval_globals):
#        time_vec = numpy.linspace(curr_time, save_time, 101)
#        occurred = []
#
#        sol_vec = scipy.integrate.odeint(BasicHybridSolver.f, y0,
#                                                    t=time_vec,
#                                                    args=(curr_state, model.listOfSpecies, model.listOfParameters, model.listOfReactions, model.listOfRateRules, eval_globals)
#                                                    )
#        check_time = len(time_vec) - 1
#        last_check_time = None
#        # see if the timepoint had any reactions over 0
#        # start from the end, look at all time points, if
#        while True:
#            for i, r in enumerate(model.listOfReactions):
#                if sol_vec[check_time,i] > 0:
#                    occurred.append(r)
#                n_occur = len(occurred)
#
#            if n_occur == 1:
#                break
#            elif n_occur > 1:
#                #decrease check_time
#                last_check_time = check_time
#                check_time -= 1
#                occurred = []
#            elif n_occur == 0:
#                #increase check_time
#                if last_check_time is not None and last_check_time == check_time + 1:
#                    # recurse down, simulate a much finer scale (1/100 the time)
#                    # UPDATE THE STATE of the continuous species
#                    for i, s in enumerate(model.listOfRateRules):
#                        curr_state[s] = sol_vec[check_time, i+len(model.listOfReactions)]
#                    return  get_reaction(curr_state, sol_vec[check_time,:], model, time_vec[check_time], time_vec[last_check_time], eval_globals)
#                last_check_time = check_time
#                check_time += 1
#                occurred = []
#            if check_time <= 0:
#                # recurse down, simulate a much finer scale (1/100 the time)
#                return  get_reaction(curr_state, y0, model, curr_time, time_vec[1], eval_globals)
#            if check_time >= len(time_vec):
#                check_time = len(time_vec)-1
#                break
#
#        # UPDATE THE STATE of the continuous species
#        for i, s in enumerate(model.listOfRateRules):
#            curr_state[s] = sol_vec[check_time, i+len(model.listOfReactions)]
#
#        print("Reaction Fired: ", occurred)
#        return occurred[0], sol_vec[check_time,:], curr_state, time_vec[check_time]


    @staticmethod
    def get_reaction(curr_state, y0, model, curr_time, save_time, eval_globals, debug):
        multiple = False  # flag variable for multiple reactions
        current = None      #current matrix state of species
        rhs = ode(BasicHybridSolver.f) #set function as ODE object
        rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model.listOfSpecies, model.listOfParameters, model.listOfReactions, model.listOfRateRules, eval_globals)
        #TODO: do we need to pass species and parameters to "f()"??
        last_state = y0
        last_time = curr_time
        step = max(0.1, save_time / 10)

        while rhs.successful():
            # Save previous state and time
            if current is not None and not multiple:    # NO reaction Fired.
                last_state = current
                last_time = int_time

            int_time = last_time+step

            if int_time > save_time:
                int_time = save_time

            if debug:
                print("Curr Time: ", curr_time, " Save time: ", save_time, " Int Time: ", int_time)
            current = rhs.integrate(int_time)   # current holds integration from current_time to int_time

            occurred = []
            for i, r in enumerate(model.listOfReactions):
                if current[i] > 0:
                    occurred.append(r)
            n_occur = len(occurred)
            if n_occur == 1:
                break
            elif n_occur > 1:
                multiple = True
                step = step * .5
                curr_time = last_time
                rhs = ode(BasicHybridSolver.f)
                rhs.set_initial_value(last_state, last_time).set_f_params(curr_state, model.listOfSpecies,
                                                                          model.listOfParameters, model.listOfReactions, model.listOfRateRules, eval_globals)
            elif int_time == save_time:
                occurred.append(None)
                break
            else:
                multiple = False


        # UPDATE THE STATE of the continuous species
        for i, s in enumerate(model.listOfRateRules):
            curr_state[s] = current[i+len(model.listOfReactions)]

        if debug:
            print("Reaction Fired: ", occurred)
            print(current)
        return occurred[0], current, curr_state, int_time

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
             **kwargs):
        if not isinstance(self, BasicHybridSolver):
            self = BasicHybridSolver()
       
        if debug:
            print("t = ",t)
            print("increment = ", increment)
        
        import math
        eval_globals = math.__dict__
        
        
        num_deterministic = 0
        #TODO
        #num_deterministic must change when toggling deterministic
        y0 = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
        propensities = {}
        curr_state = {}
        curr_time = 0
        curr_state['vol'] = model.volume
        save_time = 0

        results = {'time': []}

        for s in model.listOfSpecies:
            # initialize populations
            curr_state[s] = model.listOfSpecies[s].initial_value
            results[s] = []

        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        propensity_sum = 0
        for i, r in enumerate(model.listOfReactions):
            y0[i] = (math.log(random.uniform(0, 1)))

        while save_time < t:
            while curr_time < save_time:
                for i, rr in enumerate(model.listOfRateRules):
                    spec = model.listOfRateRules[rr].species.name
                    y0[i + len(model.listOfReactions)] = curr_state[spec]

                for i, r in enumerate(model.listOfReactions):
                    propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                    if debug:
                        print("Propensity of ", r, " is ", propensities[r])
                    propensity_sum += propensities[r]
                if propensity_sum <= 0:
                    while save_time <= t:
                        results['time'].append(save_time)
                        for s in model.listOfSpecies:
                            results[s].append(curr_state[s])
                        save_time += increment
                    return results

                #TODO: change ".1" to Salis et al. eq (16)
                reaction, y0, curr_state, curr_time = self.get_reaction(curr_state, y0, model, curr_time, save_time, eval_globals, debug=debug)

                # Update curr_state with the result of the SSA reaction that fired
                if reaction is not None:
                    for i, r in enumerate(model.listOfReactions):
                        if r == reaction:
                            y0[i] = (math.log(random.uniform(0, 1)))
                            break
                    for reactant in model.listOfReactions[reaction].reactants:
                        curr_state[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]
                    for product in model.listOfReactions[reaction].products:
                        curr_state[str(product)] += model.listOfReactions[reaction].products[product]
            results['time'].append(save_time)
            for i, s in enumerate(model.listOfSpecies):
                results[s].append(curr_state[s])
            save_time += increment
        return results
