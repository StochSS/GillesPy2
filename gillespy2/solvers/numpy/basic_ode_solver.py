from gillespy2.core import GillesPySolver
from scipy.integrate import odeint
import numpy as np


class BasicODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "BasicODESolver"
    @staticmethod
    def rhs(y0, t, species, parameters, reactions):
        """
        The right hand side of the differential equation, uses scipy.integrate odeint
        :param y0: state as a list
        :param t: time as a numpy array
        :param species: model list of species
        :param parameters: model list of parameters
        :param reactions: model list of reactions
        :return: integration step
        """
        curr_state = {}
        state_change = {}
        #   TODO    Fix Volume to take input volume
        curr_state['vol'] = 1
        propensity = {}
        results = []

        for i, s in enumerate(species):
            curr_state[s] = y0[i]
            state_change[s] = 0

        for p in parameters:
            curr_state[p] = parameters[p].value

        for r in reactions:
            propensity[r] = eval(reactions[r].propensity_function, curr_state)  # assumption that prop is massAction
            for react in reactions[r].reactants:
                state_change[str(react)] -= propensity[r]
            for prod in reactions[r].products:
                state_change[str(prod)] += propensity[r]

        for s in species:
            results.append(state_change[s])

        return (results)

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, profile=False, show_labels=True, **kwargs):
        """

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.  This is deterministic and will always have same results
        :param increment: time step increment for plotting
        :param seed: random seed, has no effect
        :param debug: not implemented
        :param profile: not implemented
        :param show_labels: not implemented
        :param kwargs:
        :return:
        """

        if show_labels:
            results = []
        else:
            num_save_times = int((t / increment)) + 1
            results = np.empty((number_of_trajectories, num_save_times, (len(model.listOfSpecies)+1)))
        for traj_num in range(number_of_trajectories):
            y0 = []
            for s in model.listOfSpecies:
                y0.append(model.listOfSpecies[s].initial_value)
            time = np.arange(0, t, increment)
            result = odeint(BasicODESolver.rhs, y0, time,
                         args=(model.listOfSpecies, model.listOfParameters, model.listOfReactions))

            #   TODO: Optimize show_labels
            if show_labels:
                results_as_dict = {}
                results_as_dict['time'] = []
                for i in range(len(time)):
                    results_as_dict['time'].append(time[i])
                for i, s in enumerate(model.listOfSpecies):
                    results_as_dict[s] = []
                    for j in range(len(result)):
                        results_as_dict[s].append(result[j][i])
                results.append(results_as_dict)
            else:
                for i in range(len(time)):
                    results[traj_num, i, 0] = time[i]
                for i, s in enumerate(model.listOfSpecies):
                    for j in range(len(result)):
                        results[traj_num, j, i+1] = result[j, i]

        return results
        # return[results, time]
        # return np.append(results, time.reshape(len(time), 1), axis=1)
