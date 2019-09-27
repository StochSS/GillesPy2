"""GillesPy2 Solver for ODE solutions."""

from scipy.integrate import ode
from scipy.integrate import odeint
from collections import OrderedDict
import numpy as np
from gillespy2.core import GillesPySolver, log


class BasicODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "BasicODESolver"
    
    def __init__(self):
        name = "BasicODESolver"

    @staticmethod
    def __f(t, y, curr_state, model, c_prop):
        """
        The right hand side of the differential equation, uses scipy.integrate odeint
        :param t: time as a numpy array
        :param y: species pops as a list
        :param current_state: dictionary of eval variables
        :param model: model being simulated
        :param c_prop: precompiled reaction propensity function
        :return: integration step
        """
        curr_state['t'] = t
        state_change =  OrderedDict()
        for i, species in enumerate(model.listOfSpecies):
            curr_state[species] = y[i]
            state_change[species] = 0
        propensity = OrderedDict()
        for r_name, reaction in model.listOfReactions.items():
            propensity[r_name] = eval(c_prop[r_name], curr_state)
            for react, stoich in reaction.reactants.items():
                state_change[str(react)] -= propensity[r_name] * stoich
            for prod, stoich in reaction.products.items():
                state_change[str(prod)] += propensity[r_name] * stoich
        state_change = list(state_change.values())
        return state_change

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, 
            show_labels=True, integrator='lsoda', integrator_options={}, **kwargs):
        """

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.
            This is deterministic and will always have same results
        :param increment: time step increment for plotting
        :param show_labels: If true, simulation returns a list of trajectories, where each list entry is a dictionary containing key value pairs of species : trajectory.  If false, returns a numpy array with shape [traj_no, time, species]
        :param integrator: integrator to be used form scipy.integrate.ode. Options include 'vode', 'zvode', 'lsoda', 'dopri5', and 'dop835'.  For more details, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        :param integrator_options: a dictionary containing options to the scipy integrator. for a list of options, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}
        :param kwargs:
        :return:
        """
        if not isinstance(self, BasicODESolver):
            self = BasicODESolver()
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. Consider running with only 1 trajectory.")

        start_state = [model.listOfSpecies[species].initial_value for species in model.listOfSpecies]

        # create mapping of species dictionary to array indices
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        # create numpy array for timeline
        timeline = np.linspace(0, t, round(t / increment + 1))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.empty((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        # copy initial populations to base
        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value

        # compile reaction propensity functions for eval
        c_prop = OrderedDict()
        for r_name, reaction in model.listOfReactions.items():
            c_prop[r_name] = compile(reaction.ode_propensity_function, '<string>', 'eval')

        result = trajectory_base[0]
        curr_time = 0
        entry_count = 0

        y0 = [0] * len(model.listOfSpecies)
        curr_state = OrderedDict()
        for i, s in enumerate(model.listOfSpecies.values()):
            curr_state[s.name] = s.initial_value
            y0[i] = s.initial_value
        for p_name, param in model.listOfParameters.items():
            curr_state[p_name] = param.value
        rhs = ode(BasicODESolver.__f).set_integrator(integrator, **integrator_options)
        rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model, c_prop)

        while entry_count < timeline.size - 1:
            int_time = curr_time + increment
            entry_count += 1
            y0 = rhs.integrate(int_time)
            curr_time += increment
            for i, spec in enumerate(model.listOfSpecies):
                curr_state[spec] = y0[i]
                result[entry_count][i+1] = curr_state[spec]

        if show_labels:
            results_as_dict = {
                'time': timeline
            }
            for i, species in enumerate(model.listOfSpecies):
                results_as_dict[species] = result[:, i+1]
            results = [results_as_dict] * number_of_trajectories
        else:
            results = np.stack([result] * number_of_trajectories, axis=0)

        return results
