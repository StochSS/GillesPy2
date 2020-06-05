"""GillesPy2 Solver for ODE solutions."""

from threading import Thread, Event
from scipy.integrate import ode
from scipy.integrate import odeint
from collections import OrderedDict
import numpy as np
from gillespy2.core import GillesPySolver, log, gillespyError


class BasicODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "BasicODESolver"
    rc = 0
    stop_event = None
    result = None
    pause_event = None
    
    def __init__(self):
        name = "BasicODESolver"
        rc = 0
        stop_event = None
        pause_event = None
        result = None

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
                state_change[react.name] -= propensity[r_name] * stoich
            for prod, stoich in reaction.products.items():
                state_change[prod.name] += propensity[r_name] * stoich
        state_change = list(state_change.values())
        return state_change

    @classmethod
    def get_solver_settings(self):
        """
        :return: Tuple of strings, denoting all keyword argument for this solvers run() method.
        """
        return ('model', 't', 'number_of_trajectories', 'increment', 'show_labels', 'integrator', 'integrator_options'
                'timeout')

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, 
            show_labels=True, integrator='lsoda', integrator_options={}, 
            timeout=None, resume = None, **kwargs):
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
        :param resume: Result of a previously run simulation, to be resumed
        :return:
        """
        if isinstance(self, type):
            self = BasicODESolver()
        self.stop_event = Event()
        self.pause_event = Event()

        if timeout is not None and timeout <=0: timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. Consider running with only 1 trajectory.")
        sim_thread = Thread(target=self.___run, args=(model,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'increment':increment, 'show_labels':show_labels, 
                                        'timeout':timeout, 'resume':resume,'integrator':integrator,
                                                                       'integrator_options':integrator_options})
        try:
            sim_thread.start()
            sim_thread.join(timeout=timeout)
            self.stop_event.set()
            while self.result is None: pass
        except:
            self.pause_event.set()
            while self.result is None: pass
        if hasattr(self, 'has_raised_exception'):
            raise self.has_raised_exception
        return self.result, self.rc

    def ___run(self, model, t=20, number_of_trajectories=1, increment=0.05, timeout=None,
            show_labels=True, integrator='lsoda', integrator_options={}, resume=None, **kwargs):

        try:
            self.__run(model, t, number_of_trajectories, increment, timeout,
                        show_labels, integrator, integrator_options, resume, **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, model, t=20, number_of_trajectories=1, increment=0.05, timeout=None,
            show_labels=True, integrator='lsoda', integrator_options={}, resume=None, **kwargs):

        timeStopped = 0
        if not (resume is None):
            if show_labels == False:
                if t < resume[0][-1][0]:
                    raise gillespyError.ExecutionError(
                        "'t' must be greater than previous simulations end time, or set in the run() method as the "
                        "simulations next end time")
            else:
                if t < resume['time'][-1]:
                    raise gillespyError.ExecutionError(
                        "'t' must be greater than previous simulations end time, or set in the run() method as the "
                        "simulations next end time")

        start_state = [model.listOfSpecies[species].initial_value for species in model.listOfSpecies]

        # create mapping of species dictionary to array indices
        species_mappings = model.sanitized_species_names()
        species = list(species_mappings.keys())
        parameter_mappings = model.sanitized_parameter_names()
        number_species = len(species)

        if not (resume is None):
            # start where we last left off if resuming a simulation
            if show_labels == False:
                lastT = resume[0][-1][0]
                step = lastT - resume[0][-2][0]
            else:
                lastT = resume['time'][-1]
                step = lastT - resume['time'][-2]
            timeline = np.arange(lastT, t+step, step)
        else:
            timeline = np.linspace(0, t, int(round(t / increment + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        # copy initial populations to base
        if not (resume is None):
            tmpSpecies = {}
            # Set initial values of species to where last left off
            if show_labels == False:
                for i, s in enumerate(species):
                    tmpSpecies[s] = resume[0][-1][i + 1]
            else:
                for i in species:
                    tmpSpecies[i] = resume[i][-1]
            for i, s in enumerate(species):
                trajectory_base[:, 0, i + 1] = tmpSpecies[s]
        else:
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

        if not (resume is None):
            for i,s in enumerate(tmpSpecies):
                curr_state[s] = tmpSpecies[s]
                y0[i] = tmpSpecies[s]
        else:
            for i, s in enumerate(model.listOfSpecies.values()):
                curr_state[s.name] = s.initial_value
                y0[i] = s.initial_value

        for p_name, param in model.listOfParameters.items():
            curr_state[p_name] = param.value
        rhs = ode(BasicODESolver.__f).set_integrator(integrator, **integrator_options)
        rhs.set_initial_value(y0, curr_time).set_f_params(curr_state, model, c_prop)

        while entry_count < timeline.size - 1:
            if self.stop_event.is_set():
                self.rc = 33
                break
            if self.pause_event.is_set():
                timeStopped = timeline[entry_count]
                break

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

        if show_labels == False and timeStopped != 0:
            cutoff = np.where(results[0][:, 0] == timeStopped)
            #Find where index is of timestopped. Ex, timestopped @50
            #index of time 50 could be 4,0, 4th row, 0'th index
            results = np.array([results[0][:int(cutoff[0])]])
        elif timeStopped != 0 and show_labels != False:
            if timeStopped != results[0]['time'][-1]:
                tester = np.where(results[0]['time'] > timeStopped)[0].size
                index = np.where(results[0]['time'] == timeStopped)[0][0]
            if tester > 0:
                for i in results[0]:
                    results[0][i] = results[0][i][:index]

        if not (resume is None):
            # If resuming, combine old pause with new data, and delete any excess null data
            if show_labels == False:
                resume = np.array([resume[0][:-1]])
                results = np.array(np.concatenate((resume, results), axis=1))
            else:
                for i in results[0]:
                    oldData = resume[i][:-1]
                    newData = results[0][i]
                    results[0][i] = np.concatenate((oldData, newData), axis=None)

        self.result = results
        return results, self.rc
