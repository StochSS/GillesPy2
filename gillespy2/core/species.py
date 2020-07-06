from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.gillespyError import *
import numpy as np

class Species(SortableObject):
    """
    Chemical species. Can be added to Model object to interact with other
    species or time.

    :param name: The name by which this species will be called in reactions and within the model.
    :type name: str
    :param initial_value: Initial population of this species. If this is not provided as an int,
    the type will be changed when it is added by numpy.int
    :type initial_value: int >= 0
    :param constant: If true, the value of the species cannot be changed (currently TauHybridSolver only)
    :type constant: bool
    :param boundary_condition: If true, species can be changed by events and rate rules, but not by reactions.
    (TauHybridSolver only)
    :type boundary_condition: bool
    :param mode: ***FOR USE WITH BasicTauHybridSolver ONLY***
    Sets the mode of representation of this species for the TauHybridSolver,
    can be discrete, continuous, or dynamic.
    mode='dynamic' - Allows a species to be represented as either discrete or continuous
    mode='continuous' - Species will only be represented as continuous
    mode='discrete' - Species will only be represented as discrete
    :type mode: str
    :param allow_negative_populations: If true, population can be reduces below 0.
    :type allow_negative_populations: bool
    :param switch_tol: ***FOR USE WITH BasicTauHybridSolver ONLY***
    Tolerance level for considering a dynamic species deterministically, value is compared to an estimated sd/mean
    population of a species after a given time step. This value will be used if a switch_min is not provided.
    The default value is 0.03
    :type switch_tol: float
    :param switch_min:  ***FOR USE WITH BasicTauHybridSolver ONLY***
    Minimum population value at which species will be represented as continuous. If a value is given, switch_min will be
    used instead of switch_tol
    :type switch_min: float
    """

    def __init__(self, name="", initial_value=0, constant=False,
                 boundary_condition=False, mode=None,
                 allow_negative_populations=False, switch_min=0,
                 switch_tol=0.03):
        # A species has a name (string) and an initial value (positive integer)
        self.name = name
        self.constant = constant
        self.boundary_condition = boundary_condition
        self.mode = mode
        self.allow_negative_populations = allow_negative_populations
        self.switch_min = switch_min
        self.switch_tol = switch_tol

        mode_list = ['continuous', 'dynamic', 'discrete', None]

        if self.mode not in mode_list:
            raise SpeciesError('Species mode must be either \'continuous\', \'dynamic\', \'discrete\', or '
                               '\'unspecified(default to dynamic for BasicTauHybridSolver)\'.')
        if mode == 'continuous':
            self.initial_value = np.float(initial_value)
        else:
            if np.int(initial_value) != initial_value:
                raise ValueError(
                    "'initial_value' for Species with mode='discrete' must be an integer value. Change to "
                    "mode='continuous' to use floating point values.")
            self.initial_value = np.int(initial_value)
        if not allow_negative_populations:
            if self.initial_value < 0:
                raise ValueError('A species initial value must be non-negative unless allow_negative_populations=True')

    def __str__(self):
        print_string = self.name
        print_string += ': ' + str(self.initial_value)
        '''
        print_string += '\n\tInitial Value: ' + str(self.initial_value)
        print_string += '\n\tConstant: ' + str(self.constant)
        print_string += '\n\tBoundary Condition: ' + str(self.boundary_condition)
        print_string += '\n\tMode: ' + self.mode
        print_string += '\n\tAllow Negative Populations: ' + str(self.allow_negative_populations)
        '''
        return print_string

    def set_initial_value(self, num):
        """
        Setter method for initial_value of a population
        :param num: Integer to set initial species population
        :raises SpeciesError: If num is non-negative or a decimal number
        """
        if isinstance(num, float) and (self.mode != 'dynamic' or self.mode != 'continuous'):
            raise SpeciesError("Mode set to discrete, species must be an integer number.")
        if num < 0 and self.allow_negative_populations == False:
            raise SpeciesError("Species population must be non-negative, or allow_negative_populations "
                               "must be set to True")
        self.initial_value = num