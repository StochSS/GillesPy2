mport gillespy2
from .gillespySolver import GillesPySolver
import random
import math
import numpy
from scipy.integrate import ode


class BasicRootSolver(GillesPySolver):

    @staticmethod
    def f(t, y, pops, species, parameters, reactions):
        curr_state = {}
        curr_state['vol'] = 1
        state_change = []
        propensity = {}
        results = []
        for i, s in enumerate(species):
            curr_state[s] = pops[i]
        for p in parameters:
            curr_state[p] = parameters[p].value
        for i, r in enumerate(reactions):
            state_change.append(eval(reactions[r].propensity_function, curr_state))
        return state_change

    def get_reaction(self, populations, y0, current_time, model, step):
        mult = False  # flag variable for multiple reactions
        int_time = 0
        current = None
        current_time = 0

        rhs = ode(BasicRootSolver.f)
        rhs.set_initial_value(y0, current_time).set_f_params(populations, model.listOfSpecies, model.listOfParameters,
                                                             model.listOfReactions)
        last_state = y0
        last_time = current_time

        while rhs.successful():

            if current is not None and not mult:
                last_state = current
                last_time = int_time

            int_time += step
            current = rhs.integrate(int_time)

            occurred = []
            for i, r in enumerate(model.listOfReactions):
                if (current[i] > 0):
                    occurred.append(r)
            n_occur = len(occurred)

            if (n_occur == 1):
                reaction = occurred[0]
                break
            elif (n_occur > 1):
                mult = mult * (-1)
                if (mult == True):
                    mult = False
                mult = True
                step = step * (.5)
                int_time = last_time
                rhs = ode(BasicRootSolver.f)
                rhs.set_initial_value(last_state, last_time).set_f_params(populations, model.listOfSpecies,
                                                                          model.listOfParameters, model.listOfReactions)
            else:
                mult = False

        return reaction

    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1, increment=0.05, seed=None, debug=False, show_labels=False,
            stochkit_home=None):
        self = BasicRootSolver()
        y0 = [0] * len(model.listOfReactions)
        populations = []
        propensities = {}
        curr_state = {}
        curr_time = 0
        rhs = ode(BasicRootSolver.f)
        curr_state['vol'] = 1
        save_time = 0

        endTime = t
        results = {'time': []}

        testCount = 0
        step = increment

        for s in model.listOfSpecies:
            populations.append(model.listOfSpecies[s].initial_value)  # initialize populations
            curr_state[s] = model.listOfSpecies[s].initial_value
            results[s] = []

        for p in model.listOfParameters:
            curr_state[p] = model.listOfParameters[p].value

        while (curr_time < endTime):
            testCount += 1
            propsum = 0
            for i, r in enumerate(model.listOfReactions):
                y0[i] = (math.log(random.uniform(0, 1)))
                propensities[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                propsum += propensities[r]
            if (propsum <= 0):
                while (save_time <= t):
                    results['time'].append(save_time)
                    for s in model.listOfSpecies:
                        results[s].append(curr_state[s])
                    save_time += increment
                return results
            step = increment


            reaction = self.get_reaction(populations, y0, curr_time, model, .1)

            for i, r in enumerate(model.listOfReactions):
                if (r == reaction):
                    break

            curr_time += -1 * y0[i] / propensities[reaction]
            while save_time < curr_time <= t:
                results['time'].append(save_time)
                for s in model.listOfSpecies:
                    results[s].append(curr_state[s])
                save_time += increment
            # apply reaction

            for reactant in model.listOfReactions[reaction].reactants:
                curr_state[str(reactant)] -= model.listOfReactions[reaction].reactants[reactant]
            for product in model.listOfReactions[reaction].products:
                curr_state[str(product)] += model.listOfReactions[reaction].products[product]
            for i, s in enumerate(model.listOfSpecies):
                # results[s].append(curr_state[s])
                populations[i] = curr_state[s]

        return results

    # else:
    # TODO: need to account for 'show_labels'
