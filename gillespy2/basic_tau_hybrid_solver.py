from .gillespySolver import GillesPySolver
import random
from scipy.integrate import ode
import math
eval_globals = math.__dict__


class BasicTauHybridSolver(GillesPySolver):
    name = "Basic Tau Hybrid Solver"

    def __init__(self, debug=False):
        self.debug = debug
    
    
    @staticmethod
    def f(t, y, curr_state, reactions, rate_rules, propensities):
        '''
        Evaluate the propensities for the reactions and the RHS of the RateRules.
        '''
        curr_state['t'] = t
        state_change = []
        for i, r in enumerate(reactions):
            propensities[r] = eval(reactions[r].propensity_function, eval_globals, curr_state)
            state_change.append(propensities[r])
        for i, rr in enumerate(rate_rules):
            state_change.append(eval(rate_rules[rr].expression,  eval_globals, curr_state))

        return state_change

    @staticmethod
    def get_reaction_integrate(step, curr_state, y0, model, curr_time, propensities):
        ''' Helper function to perform the ODE integration of one step '''
        rhs = ode(BasicTauHybridSolver.f)  # set function as ODE object
        rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model.listOfReactions,
                                                          model.listOfRateRules, propensities)
        # rhs.set_integrator('dop853')

        current = rhs.integrate(step+curr_time)   # current holds integration from current_time to int_time
        #print("time started: ", curr_time, " step taken: ", step, " current: ", current)
        if rhs.successful():
            return current, curr_time + step
        else:
            # if step is < 1e-15, take a Forward-Euler step for all species ('propensites' and RateRules)
            print("NOT RHS NOT SUCCESSFUL")
            # TODO The RateRule linked species should still contain the correct value in current, verify this
            for i, rr in enumerate(model.listOfRateRules):
                print("RHS FAILED: value of continuous species is: ", current[i+len(model.listOfReactions)])
            exit(0)

    def get_reactions(self, step, curr_state, y0, model, curr_time, save_time,
                      propensities, debug):
        ''' Get the time to the next reaction by integrating the SSA reaction functions
            along with the RateRules.  Update population of species governed by rate rules
        '''

        last_state = y0
        last_time = curr_time

        if curr_time+step > save_time:
            if debug:
                print("Step exceeds save_time, changing step size from ", step,
                      " to ", save_time - curr_time)
            step = save_time - curr_time

        if debug:
            print("Curr Time: ", curr_time, " Save time: ", save_time,  "step: ", step)

        current, curr_time = self.get_reaction_integrate(step,  curr_state, y0, model,
                                                         curr_time, propensities)

        rxn_count = {}
        fired = False
        for i, r in enumerate(model.listOfReactions):
            #urn = (math.log(random.uniform(0, 1)))
            rxn_count[r] = 0
            rj = current[i]
            #print("current[i]: ", current[i], " y0[i]: ", y0[i], " rj: ", rj)
            #print(r, " rj is ", rj)
            while rj > 0:
                #print(r, " fired")
                if not fired:
                    fired = True
                rxn_count[r] += 1
                urn = (math.log(random.uniform(0, 1)))
                rj += urn

        # occurred = []
        # for i, r in enumerate(model.listOfReactions):
        #     if current[i] >= 0:
        #         occurred.append(r)
        # n_occur = len(occurred)
        # if n_occur == 1:
        #     break
        # elif n_occur > 1:
        #     if debug:
        #         print("Multiple reactions fired in this step (n=", n_occur, ") changing step size from ", step,
        #               " to ", step*0.75, "recursion_counter: ", recursion_counter)
        #     # reset state, and try again
        #     step = step * .75
        #     curr_time = last_time
        #     y0 = last_state
        #     recursion_counter += 1
        #     if recursion_counter > 20:
        #         raise Exception("get_reaction() failed, too many step size reductions: halved {0} times"
        #                         .format(recursion_counter))
        # elif curr_time >= save_time:
        #     occurred.append(None)
        #     break
        # else:
        #     # ODE was successful, but no reactions fired, advance time
        #     last_state = current
        #     last_time = curr_time
        if not fired:
            if debug:
                print("No reactions fired in this step changing step size from ", step,
                      " to ", step * 1.25)
            step = step * 1.25

        # UPDATE THE STATE of the continuous species
        for i, s in enumerate(model.listOfRateRules):
            curr_state[s] = current[i+len(model.listOfReactions)]

        if debug:
            print("Reactions Fired: ", rxn_count)
            print("y(t) = ", current)

        return rxn_count, current, curr_state, curr_time

    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
            **kwargs):
        if not isinstance(self, BasicTauHybridSolver):
            self = BasicTauHybridSolver()
        if debug:
            print("t = ", t)
            print("increment = ", increment)

        random.seed(seed)

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
        for i, r in enumerate(model.listOfReactions):   # set reactions to uniform random number and add to y0
            y0[i] = (math.log(random.uniform(0, 1)))
            if debug:
                print("Setting Random number ", y0[i], " for ", model.listOfReactions[r].name)

        while save_time < t:
            while curr_time < save_time:
                projected_reaction = None
                tau_step = None
                tau_j = {}
                for i, rr in enumerate(model.listOfRateRules):  # Add continuous species to y0
                    spec = model.listOfRateRules[rr].species.name
                    y0[i + len(model.listOfReactions)] = curr_state[spec]
                for i, r in enumerate(model.listOfReactions):
                    propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                    propensity_sum += propensities[r]
                    # Salis et al. eq (16)
                    # TODO: this needs to be optimized.  Going too big is expensive, too small is also expensive
                    if propensities[r] > 0:
                        tau_j[r] = -y0[i] / propensities[r]
                        if debug:
                            print("Propensity of ", r, " is ", propensities[r], "tau_j is ", tau_j[r])
                        if tau_step is None or tau_j[r] < tau_step:
                            tau_step = max(tau_j[r], 1e-13)
                            projected_reaction = model.listOfReactions[r]
                    else:
                        if debug:
                            print("Propensity of ", r, " is 0")
                if tau_step is None:
                    tau_step = save_time - curr_time
                if debug:
                    if projected_reaction is None:
                        print("NO projected reaction")
                    else:
                        print("Projected reaction is: ", projected_reaction.name, " at time: ", curr_time+tau_step,
                          " step size", tau_step)

                reactions, y0, curr_state, curr_time = self.get_reactions(
                    tau_step, curr_state, y0, model, curr_time, save_time, propensities, debug)
                # Update curr_state with the result of the SSA reaction that fired
                for i, r in enumerate(model.listOfReactions):
                    #print("at index: ", i, " checking ", r, ": ", reactions[r])
                    if reactions[r] > 0:
                        #print(r, " is greater than 0")
                        for reactant in model.listOfReactions[r].reactants:
                            #print("Updating reactant: ", reactant)
                            for j in range(reactions[r]):
                                #print("decrementing ", reactant)
                                curr_state[str(reactant)] -= model.listOfReactions[r].reactants[reactant]
                                #print("curr state of ", reactant, " is ", curr_state[str(reactant)])
                        for product in model.listOfReactions[r].products:
                            #print("Updating product: ", product)
                            for j in range(reactions[r]):
                                #print("incrementing", product)
                                curr_state[str(product)] += model.listOfReactions[r].products[product]
                                #print("curr state of ", product, " is ", curr_state[str(product)])
                        y0[i] += (math.log(random.uniform(0, 1)))
                        if debug:
                            print("Setting Random number ", y0[i], " for ", model.listOfReactions[r].name)
                        break

                # for reactant in model.listOfReactions[reaction].reactants:
                #     curr_state[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]
                # for product in model.listOfReactions[reaction].products:
                #         curr_state[str(product)] += model.listOfReactions[reaction].products[product]
            results['time'].append(save_time)
            for i, s in enumerate(model.listOfSpecies):
                results[s].append(curr_state[s])
            save_time += increment
        return results
