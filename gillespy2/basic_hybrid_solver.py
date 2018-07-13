from .gillespySolver import GillesPySolver
import random
import math
from scipy.integrate import ode


class BasicHybridSolver(GillesPySolver):
    name = "Basic Hybrid Solver"

    @staticmethod
    # Create the matrix equation for integration
    def f(t, y, pops, species, parameters, reactions, rate_rules):
        curr_state = {'vol': 1}
        curr_state['t'] = t
        state_change = []
        for i, s in enumerate(species):
            curr_state[s] = pops[i]
        for p in parameters:
            curr_state[p] = parameters[p].value
        for i, r in enumerate(reactions):
            state_change.append(eval(reactions[r].propensity_function, curr_state))
        for i, rr in enumerate(rate_rules):
            state_change.append(eval(rate_rules[rr].expression, curr_state))

        return state_change

    @staticmethod
    def get_reaction(populations, y0, model, step, curr_time, save_time):
        multiple = False  # flag variable for multiple reactions
        current_time = curr_time    #integration start time
        int_time = current_time    #integration end time
        current = None      #current matrix state of species
        reaction = []       #list beginning with reaction fired, and followed by all species population states at reaction time
        print("Save time at beginning of get reaction: ", save_time)
        rhs = ode(BasicHybridSolver.f) #set function as ODE object
        rhs.set_initial_value(y0, current_time).set_f_params(populations, model.listOfSpecies, model.listOfParameters, model.listOfReactions, model.listOfRateRules)
        last_state = y0
        last_time = current_time

        while rhs.successful():
            # Save previous state and time
            if current is not None and not multiple:
                last_state = current
                last_time = int_time

            int_time += step
            print("Int Time: ", int_time)
            if int_time > save_time:
                int_time = save_time

            current = rhs.integrate(int_time) # current holds integration from current_time to int_time

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
                int_time = last_time
                rhs = ode(BasicHybridSolver.f)
                rhs.set_initial_value(last_state, last_time).set_f_params(populations, model.listOfSpecies,
                                                                          model.listOfParameters, model.listOfReactions, model.listOfRateRules)
            elif int_time == save_time:
                occurred.append(None)

                break
            else:
                multiple = False


        # APPEND THE STATE
        for i, s in enumerate(model.listOfRateRules):
            populations[s] = current[i+len(model.listOfReactions)]

        print("Reaction Fired: ", reaction)
        print(current)
        return occurred[0], current, populations, int_time

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
            stochkit_home=None, **kwargs):
        if not isinstance(self, BasicHybridSolver):
            self = BasicHybridSolver()
        num_deterministic = 0
        #TODO
        #num_deterministic must change when toggling deterministic
        y0 = [0] * (len(model.listOfReactions) + len(model.listOfRateRules))
        propensities = {}
        curr_state = {}
        curr_time = 0
        curr_state['vol'] = 1
        save_time = 0

        end_time = t
        results = {'time': []}

        test_count = 0

        for s in model.listOfSpecies:
            # initialize populations
            curr_state[s] = model.listOfSpecies[s].initial_value
            results[s] = []

        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        while curr_time < end_time:
            test_count += 1
            propensity_sum = 0
            for i, r in enumerate(model.listOfReactions):
                y0[i] = (math.log(random.uniform(0, 1)))
                propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                print("Propensity of ", r, " is ", propensities[r])
                propensity_sum += propensities[r]
            for i, rr in enumerate(model.listOfRateRules):
                spec = model.listOfRateRules[rr].species.name
                y0[i+len(model.listOfReactions)] = curr_state[spec]
            if propensity_sum <= 0:
                while save_time <= t:
                    results['time'].append(save_time)
                    for s in model.listOfSpecies:
                        results[s].append(curr_state[s])
                    save_time += increment
                return results

            reaction, curr_state, populations, curr_time = self.get_reaction(populations, y0, model, .1, curr_time, save_time)

            for i, r in enumerate(model.listOfReactions):
                if r == reaction:
                    break
            print("Save time: ", save_time, " Curr Time: ", curr_time)
            while save_time-1 < curr_time <= t:
                results['time'].append(save_time)
                for i, s in enumerate(model.listOfSpecies):
                    if model.listOfSpecies[s].deterministic:
                        curr_state[s] = reaction[i + 1]     # reaction is ordered [rx fired, listOfSpecies, time]
                    results[s].append(curr_state[s])
                save_time += increment

            if reaction is not None:
                for reactant in model.listOfReactions[reaction[0]].reactants:
                    curr_state[str(reactant)] -= model.listOfReactions[reaction[0]].reactants[reactant]
                for product in model.listOfReactions[reaction[0]].products:
                    curr_state[str(product)] += model.listOfReactions[reaction[0]].products[product]

        return results
