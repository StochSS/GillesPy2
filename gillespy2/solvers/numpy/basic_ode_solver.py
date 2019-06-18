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
        curr_state = {}
        state_change = {}
        curr_state['vol'] = model.volume

        for i, species in enumerate(model.listOfSpecies):
            curr_state[species] = start_state[i]
            state_change[species] = 0

        for parameter in model.listOfParameters:
            curr_state[parameter] = model.listOfParameters[parameter].value

        propensity = {}
        for reaction in model.listOfReactions:
            propensity[reaction] = eval(
                model.listOfReactions[reaction].propensity_function, curr_state)
            # assumption that prop is massAction
            for react in model.listOfReactions[reaction].reactants:
                state_change[str(react)] -= propensity[reaction]
            for prod in model.listOfReactions[reaction].products:
                state_change[str(prod)] += propensity[reaction]

        results = [state_change[species] for species in model.listOfSpecies]

        return results

    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, debug=False, profile=False, show_labels=True, max_steps=0, **kwargs):
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
        :param max_steps: Defaults to 0 for odeint
            When using deterministic methods, specifies the maximum number of steps permitted for each integration point in t.
        :param kwargs:
        :return:
        """
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. Consider running with only 1 trajectory.")

        start_state = [model.listOfSpecies[species].initial_value for species in model.listOfSpecies]
        timeline = np.linspace(0, t, (t // increment + 1))
        result = odeint(BasicODESolver.rhs, start_state, timeline, args=(model,), mxstep=max_steps)
        result = np.hstack((np.expand_dims(timeline, -1), result))

        if show_labels:
            results_as_dict = {
                'time': timeline
            }
            for i, species in enumerate(model.listOfSpecies):
                results_as_dict[species] = result[:, i]
            results = [results_as_dict] * number_of_trajectories
        else:
            results = np.stack([result] * number_of_trajectories, axis=0)

        return results
