# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from gillespy2.core.jsonify import Jsonify
from gillespy2.core.sortableobject import SortableObject

from gillespy2.core.gillespyError import SpeciesError

class Species(SortableObject, Jsonify):
    """
    Chemical species. Can be added to Model object to interact with other
    species or time.

    :param name: The name by which this species will be called in reactions and within the model.
    :type name: str

    :param initial_value: Initial population (discrete) or concentration (continuous) of this species.
    :type initial_value: int | float

    :param constant: If true, the value of the species cannot be changed (currently TauHybridSolver only)
    :type constant: bool

    :param boundary_condition: If true, species can be changed by events and rate rules, but not by reactions.
        (TauHybridSolver only)
    :type boundary_condition: bool

    :param mode: ***FOR USE WITH TauHybridSolver ONLY***
        Sets the mode of representation of this species for the TauHybridSolver,
        can be discrete, continuous, or dynamic.
        mode='dynamic' - Allows a species to be represented as either discrete or continuous
        mode='continuous' - Species will only be represented as continuous
        mode='discrete' - Species will only be represented as discrete
    :type mode: str

    :param allow_negative_populations: If true, population can be reduced below 0.
    :type allow_negative_populations: bool

    :param switch_tol: ***FOR USE WITH TauHybridSolver ONLY***
        Tolerance level for considering a dynamic species deterministically, value is compared to an estimated sd/mean
        population of a species after a given time step. This value will be used if a switch_min is not provided.
        The default value is 0.03
    :type switch_tol: float

    :param switch_min:  ***FOR USE WITH TauHybridSolver ONLY***
        Minimum population value at which species will be represented as continuous.
        If a value is given, switch_min will be used instead of switch_tol
    :type switch_min: float

    :raises SpeciesError: Arg is of invalid type.  Required arg set to None.  Arg value is outside of accepted bounds.
    """

    def __init__(self, name=None, initial_value=0, constant=False, boundary_condition=False, mode=None,
                 allow_negative_populations=False, switch_min=0, switch_tol=0.03):
        # A species has a name (string) and an initial value (positive integer)
        self.name = name
        self.constant = constant
        self.boundary_condition = boundary_condition
        self.mode = mode
        self.allow_negative_populations = allow_negative_populations
        self.switch_min = switch_min
        self.switch_tol = switch_tol

        if initial_value is None:
            raise SpeciesError("initial_value can't be None type.")
        if isinstance(initial_value, str):
            try:
                initial_value = float(initial_value)
            except ValueError:
                pass
        self.validate(initial_value=initial_value)

        self.initial_value = float(initial_value) if self.mode == "continuous" else int(initial_value)

    def __str__(self):
        print_string = self.name
        print_string += ': ' + str(self.initial_value)
        return print_string

    def set_initial_value(self, initial_value):
        """
        Setter method for initial_value of a population

        :param initial_value: Initial population (discrete) or concentration (continuous) of this species.
        :type initial_value: int | float

        :raises SpeciesError: initial_value is of invalid type.  initial_value set to None.  \
                              initial_value is a float when mode != 'continuous'. \
                              initial_value is negative when allow_negative_populations=False.
        """
        if initial_value is None:
            raise SpeciesError("initial_value can't be None type.")
        if isinstance(initial_value, str):
            try:
                initial_value = float(initial_value)
            except ValueError:
                pass
        self.validate(initial_value=initial_value, coverage="initial_value")

        self.initial_value = float(initial_value) if self.mode == "continuous" else int(initial_value)

    def validate(self, initial_value=None, coverage="all"):
        """
        Validate the species.

        :param initial_value: Initial population (discrete) or concentration (continuous) of this species.
        :type initial_value: int | float

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises SpeciesError: Attribute is of invalid type.  Required attribute set to None.  \
                              Attribute is value outside of accepted bounds.
        """
        # Check name
        if coverage in ("all", "name"):
            if self.name is None:
                raise SpeciesError("name can't be None type.")
            if not isinstance(self.name, str):
                raise SpeciesError(f"name must be of type str not {type(self.name)}.")
            if self.name == "":
                raise SpeciesError("name can't be an empty string.")

        # Check initial_value
        if coverage in ("all", "initial_value"):
            if initial_value is None:
                initial_value = self.initial_value

            if initial_value is None:
                raise SpeciesError("initial_value can't be None type.")
            if not isinstance(initial_value, (float, int)):
                raise SpeciesError(f"initial_value must be of type float or int not {type(initial_value)}.")
            if self.mode != "continuous" and int(initial_value) != initial_value:
                raise SpeciesError(
                    """
                    initial_value with mode='discrete' must be an integer value.
                    Change to mode='continuous' to use floating point values.
                    """
                )
            if not self.allow_negative_populations and initial_value < 0:
                raise SpeciesError(
                    'A species initial value must be non-negative unless allow_negative_populations=True'
                )

        # Check constant
        if coverage in ("all", "constant"):
            if not isinstance(self.constant, bool):
                errmsg = f"constant must be of type bool not {type(self.constant)}."
                raise SpeciesError(errmsg)

        # Check boundary_condition
        if coverage in ("all", "boundary_condition"):
            if not isinstance(self.boundary_condition, bool):
                errmsg = f"boundary_condition must be of type bool not {type(self.boundary_condition)}."
                raise SpeciesError(errmsg)

        # Check mode
        if coverage in ("all", "mode"):
            mode_list = ['continuous', 'dynamic', 'discrete', None]

            if self.mode not in mode_list:
                raise SpeciesError(
                    f"""
                    mode must be 'continuous', 'dynamic', 'discrete', or
                    unspecified (defaults to 'dynamic' for TauHybridSolver) not {self.mode}.
                    """
                )

        # Check allow_negative_populations
        if coverage in ("all", "allow_negative_populations"):
            if not isinstance(self.allow_negative_populations, bool):
                errmsg = f"allow_negative_populations must be of type bool not {type(self.allow_negative_populations)}."
                raise SpeciesError(errmsg)

        # Check switch_tol
        if coverage in ("all", "switch_tol"):
            if self.switch_tol is None:
                raise SpeciesError("switch_tol can't be None type.")
            if not isinstance(self.switch_tol, (float, int)):
                raise SpeciesError(f"switch_tol must of type float or int not {type(self.switch_tol)}")
            if self.switch_tol < 0:
                raise SpeciesError(f"switch_tol must be a positive value not {self.switch_tol}")

        # Check switch_min
        if coverage in ("all", "switch_min"):
            if self.switch_min is None:
                raise SpeciesError("switch_min can't be None type.")
            if not isinstance(self.switch_min, (float, int)):
                raise SpeciesError(f"switch_min must of type float or int not {type(self.switch_min)}")
            if self.switch_min < 0:
                raise SpeciesError(f"switch_min must be a positive value not {self.switch_min}")
