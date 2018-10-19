import gillespy2
from .gillespySolver import GillesPySolver
import random
import math


class BasicSSASolver(GillesPySolver):
    """ TODO
    """
    name = "BasicSSASolver"

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, profile=False, debug=False, show_labels=False,stochkit_home=None):

        self.simulation_data = []

        curr_state = {}
        propensity = {}
        results = {}
        steps_taken = []

        for trajectories in range(number_of_trajectories):
            # Initialize Species population
            for s in model.listOfSpecies:
                curr_state[s] = model.listOfSpecies[s].initial_value
                results[s] = []
            curr_state['vol'] = model.volume
            results['time'] = []
            current_time = 0
            save_time = 0

            for p in model.listOfParameters:
                curr_state[p] = model.listOfParameters[p].value

            while current_time < t:
                prop_sum = 0
                cumulative_sum = 0
                reaction = None
                reaction_num = None
                for r in model.listOfReactions:
                    propensity[r] = eval(model.listOfReactions[r].propensity_function, curr_state)
                    prop_sum += propensity[r]
                reaction_num = random.uniform(0, prop_sum)
                for r in model.listOfReactions:
                    cumulative_sum += propensity[r]
                    if cumulative_sum >= reaction_num:
                        reaction = r
                        break
                if prop_sum <= 0:
                    while save_time <= t:
                        results['time'].append(save_time)
                        for s in model.listOfSpecies:
                            results[s].append(curr_state[s])
                        save_time += increment
                    return results

                tau = -1*math.log(random.random())/prop_sum
                current_time += tau
                if profile:
                    steps_taken.append(tau)

                while(current_time > save_time and current_time <= t):
                    results['time'].append(save_time)
                    for s in model.listOfSpecies:
                        results[s].append(curr_state[s])
                    save_time += increment

                for react in model.listOfReactions[reaction].reactants:
                    curr_state[str(react)] -= model.listOfReactions[reaction].reactants[react]
                for prod in model.listOfReactions[reaction].products:
                    curr_state[str(prod)] += model.listOfReactions[reaction].products[prod]
            if profile:
                print(steps_taken)
                print("Total Steps Taken", len(steps_taken))

            return results

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if show_labels:
            return self.simulation_data
    # else:
    # TODO: need to account for 'show_labels'
