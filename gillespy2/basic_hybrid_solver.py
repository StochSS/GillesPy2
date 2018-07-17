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
    def f(t, y, curr_state, reactions, rate_rules, eval_globals, propensities, curr_time):
        ''' Evaluate the propensities for the reactions and the RHS of the RateRules.

        '''
        curr_state['t'] = t
        state_change = []
        for i, r in enumerate(reactions):
            propensities[r] = eval(reactions[r].propensity_function, eval_globals, curr_state)
            state_change.append(propensities[r])
        for i, rr in enumerate(rate_rules):
            state_change.append(eval(rate_rules[rr].expression,  eval_globals, curr_state))

        return state_change

    def get_reaction_integrate(self, step, curr_state, y0, model, curr_time, eval_globals, propensities):
        ''' Helper function to perform the ODE integration of one step '''
        rhs = ode(BasicHybridSolver.f) #set function as ODE object
        rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model.listOfReactions, model.listOfRateRules, eval_globals, propensities, curr_time)
        #rhs.set_integrator('dop853')

        current = rhs.integrate(step+curr_time)   # current holds integration from current_time to int_time
        if rhs.successful():
            return current, curr_time + step
        else:
            #TODO: if step is < 1e-15, take a Forward-Euler step for all species ('propensites' and RateRules)
            raise Exception("get_reaction_integrate() failed.")



    def get_reaction(self, curr_state, y0, model, step, curr_time, save_time, eval_globals, propensities, debug):
        ''' Get the time to the next reaction by integrating the SSA reaction functions
            along with the RateRules.  Update population of species governed by rate rules
        '''
        multiple = False  # flag variable for multiple reactions
        current = None      #current matrix state of species
        last_state = y0
        last_time = curr_time
        occurred = []
        recursion_counter = 0
        time_advance_flag = False

        while True: 
            if curr_time+step > save_time:
                step = save_time - curr_time

            if debug:
                print("Curr Time: ", curr_time, " Save time: ", save_time,  "step: ", step)
                
            current, curr_time = self.get_reaction_integrate(step,curr_state, y0, model, curr_time, eval_globals, propensities)

            occurred = []
            for i, r in enumerate(model.listOfReactions):
                if current[i] > 0:
                    occurred.append(r)
            n_occur = len(occurred)
            if n_occur == 1:
                break
            elif n_occur > 1:
                if debug:
                    print("Multiple reactions fired in this step (n=",n_occur,") changing stepsize from ",step," to ",step*0.5, "recursion_counter: ",recursion_counter)
                # reset state, and try again
                step = step * .5
                curr_time = last_time
                y0 = last_state
                recursion_counter += 1
                if(recursion_counter > 20):
                    raise Exception("get_reaction() failed, too many step size reductions: halved {0} times".format(recursion_counter));
            elif curr_time >= save_time:
                occurred.append(None)
                break
            else:
                # ODE was successful, but no reactions fired, advance time
                last_state = current
                last_time = curr_time
                #TODO: if you hit this block of code more than once, double the step size
                if time_advance_flag:
                    step = step * 2
                else:
                    time_advance_flag = True


        # UPDATE THE STATE of the continuous species
        for i, s in enumerate(model.listOfRateRules):
            curr_state[s] = current[i+len(model.listOfReactions)]

        if debug:
            print("Reaction Fired: ", occurred)
            print("y(t) = ",current)

        if(len(occurred) == 0):
            #TODO make a custom exception
            raise Exception("Error scipy.integrate.ode did not succeed, return code{0}".format(rhs.get_return_code()))

        return occurred[0], current, curr_state, curr_time

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

        random.seed(seed)
        
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
                #Salis et al. eq (16)
                #TODO: this needs to be optimized.  Going too big is expensive, too small is also expensive
                tau_step = None
                for i, r in enumerate(model.listOfReactions):
                    if(propensities[r] > 0):
                        tau_j = -y0[i] / propensities[r] 
                        if debug:
                            print("Propensity of ", r, " is ", propensities[r], "tau_j is ",tau_j)
                        if(tau_step is None or tau_j < tau_step):
                            tau_step = tau_j
                if tau_step is None:
                    tau_step = save_time - curr_time


                reaction, y0, curr_state, curr_time = self.get_reaction(curr_state, y0, model, tau_step, curr_time, save_time, eval_globals, propensities, debug=debug)

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
