"""GillesPy2 Solver for ODE solutions."""

from scipy.integrate import odeint
import numpy as np
from gillespy2.core import GillesPySolver
from gillespy2.core import log


class BasicODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "BasicODESolver"

    @staticmethod
    def rhs(start_state, time, model):
        """
        The right hand side of the differential equation, uses scipy.integrate odeint
        :param start_state: state as a list
        :param t: time as a numpy array
        :param model: model being simulated
        :return: integration step
        """
        #   pylint: disable=W0613, W0123
        curr_state = {}
        state_change = {}
        curr_state['vol'] = model.volume

        propensity = {}
        results = []

        for i, species in enumerate(model.listOfSpecies):
            curr_state[species] = start_state[i]
            state_change[species] = 0

        for parameter in model.listOfParameters:
            curr_state[parameter] = model.listOfParameters[parameter].value

        for reaction in model.listOfReactions:
            propensity[reaction] = eval(
                model.listOfReactions[reaction].propensity_function, curr_state)
            # assumption that prop is massAction
            for react in model.listOfReactions[reaction].reactants:
                state_change[str(react)] -= propensity[reaction]
            for prod in model.listOfReactions[reaction].products:
                state_change[str(prod)] += propensity[reaction]

        for species in model.listOfSpecies:
            results.append(state_change[species])

        return results

    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, profile=False, show_labels=True, **kwargs):
        """

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.
            This is deterministic and will always have same results
        :param increment: time step increment for plotting
        :param seed: random seed, has no effect
        :param debug: not implemented
        :param profile: not implemented
        :param show_labels: not implemented
        :param kwargs:
        :return:
        """
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. Consider running with only 1 trajectory.")
        #   pylint: disable=R0913, R0914
        if show_labels:
            results = []
        else:
            num_save_times = int((t / increment))
            results = np.empty((number_of_trajectories,
                                num_save_times, (len(model.listOfSpecies) + 1)))

        start_state = []
        for species in model.listOfSpecies:
            start_state.append(model.listOfSpecies[species].initial_value)
        time = np.arange(0., t, increment, dtype=np.float64)
        result = odeint(BasicODESolver.rhs, start_state, time, args=(model,))

        for traj_num in range(number_of_trajectories):
            if show_labels:
                results_as_dict = {
                    'time': []
                }
                for i, timestamp in enumerate(time):
                    results_as_dict['time'].append(timestamp)
                for i, species in enumerate(model.listOfSpecies):
                    results_as_dict[species] = []
                    for row in result:
                        results_as_dict[species].append(row[i])
                results.append(results_as_dict)
            else:
                for i, timestamp in enumerate(time):
                    results[traj_num, i, 0] = timestamp
                for i in enumerate(model.listOfSpecies):
                    for j in range(len(result)):
                        results[traj_num, j, i[0] + 1] = result[j, i[0]]

        return results
