"""GillesPy2 Solver for ODE solutions."""

from threading import Thread, Event
from types import FunctionType
from scipy.integrate import ode
from collections import OrderedDict
import numpy as np
from gillespy2.core import GillesPySolver, log


class BasicODESolver(GillesPySolver):
    """
    This Solver produces the deterministic continuous solution via ODE.
    """
    name = "BasicODESolver"
    rc = 0
    stop_event = None
    result = None
    
    def __init__(self):
        name = "BasicODESolver"
        rc = 0
        stop_event = None
        result = None

    @staticmethod
    def __f(t, y, propensity, sm, propensities):
        """
        The right hand side of the differential equation, uses scipy.integrate odeint
        :param t: time as a numpy array
        :param y: species pops as a list
        :param current_state: dictionary of eval variables
        :param model: model being simulated
        :param c_prop: precompiled reaction propensity function
        :return: integration step
        """
        for i in range(len(propensity)): propensities[i] = propensity[i](y)
        state_change = np.multiply(sm, propensities).sum(axis=0)
        return state_change

    @classmethod
    def run(self, model, t=20, number_of_trajectories=1, increment=0.05, 
            show_labels=True, integrator='lsoda', integrator_options={}, 
            timeout=None, **kwargs):
        """

        :param model: gillespy2.model class object
        :param t: end time of simulation
        :param number_of_trajectories: Should be 1.
            This is deterministic and will always have same results
        :param dt: time step increment for plotting
        :param show_labels: If true, simulation returns a list of trajectories, where each list entry is a dictionary containing key value pairs of species : trajectory.  If false, returns a numpy array with shape [traj_no, time, species]
        :param integrator: integrator to be used form scipy.integrate.ode. Options include 'vode', 'zvode', 'lsoda', 'dopri5', and 'dop835'.  For more details, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
        :param integrator_options: a dictionary containing options to the scipy integrator. for a list of options, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html.
            Example use: {max_step : 0, rtol : .01}
        :param kwargs:
        :return:
        """
        if isinstance(self, type):
            self = BasicODESolver()
        self.stop_event = Event()
        if timeout is not None and timeout <=0: timeout = None
        if len(kwargs) > 0:
            for key in kwargs:
                log.warning('Unsupported keyword argument to {0} solver: {1}'.format(self.name, key))
        if number_of_trajectories > 1:
            log.warning("Generating duplicate trajectories for model with ODE Solver. Consider running with only 1 trajectory.")
        sim_thread = Thread(target=self.___run, args=(model,), kwargs={'t':t,
                                        'number_of_trajectories':number_of_trajectories,
                                        'dt':increment, 'show_labels':show_labels, 
                                        'timeout':timeout,
                                        'integrator':integrator,
                                        'integrator_options':integrator_options})
        try:
            sim_thread.start()
            sim_thread.join(timeout=timeout)
            self.stop_event.set()
            while self.result is None: pass
        except:
            pass
        if hasattr(self, 'has_raised_exception'):
            raise self.has_raised_exception
        return self.result, self.rc

    def ___run(self, model, t=20, number_of_trajectories=1, dt=0.05, timeout=None,
            show_labels=True, integrator='lsoda', integrator_options={}, **kwargs):
        try:
            self.__run(model, t, number_of_trajectories, dt, timeout,
                        show_labels, integrator, integrator_options, **kwargs)
        except Exception as e:
            self.has_raised_exception = e
            self.result = []
            return [], -1

    def __run(self, model, t=20, number_of_trajectories=1, dt=0.05, timeout=None,
            show_labels=True, integrator='lsoda', integrator_options={}, **kwargs):

        # create mappings to array indices
        species_mappings = model._listOfSpecies
        species = list(species_mappings.keys())
        parameter_mappings = model._listOfParameters
        parameters = list(parameter_mappings.keys())
        number_species = len(species)
        number_parameters = len(parameters)
        reactions = model._listOfReactions
        number_reactions = len(reactions)

        vr = [''] * number_reactions
        propensity = [0] * number_reactions
        sm = np.zeros((number_reactions, number_species+number_parameters))
        for i, r in enumerate(reactions.values()):
            vr[i] = r.propensity_function.replace('vol', str(model.volume))
            for j, s in reversed(list(enumerate(species_mappings.values()))):
                if s in r.reactants:
                    sm[i, j] -= r.reactants[s]
                if s in r.products:
                    sm[i, j] += r.products[s]
                if s in vr[i]:
                    vr[i] = vr[i].replace(s, 'y[{}]'.format(j))
            for j, p in reversed(list(enumerate(parameter_mappings.values()))):
                vr[i] = vr[i].replace(p, 'y[{}]'.format(j+number_species))
        for i in range(number_reactions):
            pcode = compile('def pf(y): return {}'.format(vr[i]), '<string>', 'exec')
            pfunc = FunctionType(pcode.co_consts[0], globals(), 'pf')
            propensity[i] = pfunc


        # create numpy array for timeline
        timeline = np.linspace(0, t, int(round(t / dt + 1)))

        # create numpy matrix to mark all state data of time and species
        trajectory_base = np.zeros((number_of_trajectories, timeline.size, number_species + 1))

        # copy time values to all trajectory row starts
        trajectory_base[:, :, 0] = timeline

        # copy initial populations to base
        for i, s in enumerate(species):
            trajectory_base[:, 0, i + 1] = model.listOfSpecies[s].initial_value

        result = trajectory_base[0]
        t0 = 0
        y0 = np.zeros(len(species)+len(parameters))

        for i, s in enumerate(species):
            y0[i] = model.listOfSpecies[s].initial_value
        for i, p in enumerate(parameters):
            y0[i+len(species)] = model.listOfParameters[p].value

        propensities = np.zeros((len(vr), 1))
        rhs = ode(BasicODESolver.__f).set_integrator(integrator, **integrator_options)
        rhs.set_initial_value(y0, t0).set_f_params(propensity, sm, propensities)
        entry_count = 0

        while entry_count < timeline.size - 1:
            if self.stop_event.is_set():
                self.rc = 33
                break
            entry_count += 1
            y = rhs.integrate(rhs.t+dt)
            result[entry_count][1:] = y[:number_species]

        if show_labels:
            results_as_dict = {
                'time': timeline
            }
            for i, species in enumerate(species):
                results_as_dict[species] = result[:, i+1]
            results = [results_as_dict] * number_of_trajectories
        else:
            results = np.stack([result] * number_of_trajectories, axis=0)
        self.result = results
        return results, self.rc
