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

import uuid

from gillespy2.core.parameter import Parameter
from gillespy2.core.species import Species
from gillespy2.core.jsonify import Jsonify

from gillespy2.core.gillespyError import EventError

class EventAssignment(Jsonify):
    """
    An EventAssignment describes a change to be performed to the current model
    simulation.  This is assignment can either be fired at the time its
    associated trigger changes from false to true, or after a specified delay,
    depending on how the Event to which it is assigned is configured.

    :param variable: Target model component to be modified by the EventAssignment
        expression. Valid target variables include gillespy2 Species,
        Parameters, and Compartments.
    :type variable: gillespy2.Species, gillespy2.Parameter

    :param expression: String to be evaluated when the event is fired.  This expression must
        be evaluable within the model namespace, and the results of it's
        evaluation will be assigned to the EventAssignment variable.
    :type expression: str

    """
    def __init__(self, name=None, variable=None, expression=None):

        if name in (None, ""):
            name = f'evn{uuid.uuid4()}'.replace('-', '_')
        else:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning("EventAssignment.name has been deprecated.")
        self.__name_deprecated = name

        self.variable = variable
        self.expression = expression

        if expression is not None:
            self.expression = str(expression)

        #TODO: ADD Compartment to valid variable types once implemented
        valid_variable_types = [Species, Parameter, str]

        if not type(variable) in valid_variable_types:
            print(variable)
            print(type(variable))
            raise EventError(
                'GillesPy2 Event Assignment variable must be a valid gillespy2 species')
        if not isinstance(self.expression, str):
            raise EventError(
                             'GillesPy2 Event Assignment expression requires a '
                             'valid string expression')

    def __str__(self):
        return f"{self.variable}: {self.expression}"

    def __getattr__(self, key):
        if key == 'name':
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning('EventAssignment.name has been deprecated.')
            return self.__name_deprecated
        raise AttributeError

class EventTrigger(Jsonify):
    """
    Trigger detects changes in model/environment conditions in order to fire an
    event.  A Trigger contains an expression, a mathematical function which can
    be evaluated to a boolean value within a model's namespace.  Upon
    transitioning from 'false' to 'true', this trigger will cause the immediate
    execution of an event's list of assignments if no delay is present, otherwise,
    the delay evaluation will be initialized.

    :param expression: String for a function calculating EventTrigger values. Should be evaluable
        in namespace of Model.
    :type expression: str

    :param value: Value of EventTrigger at simulation start, with time t=0
    :type value: bool

    :param persistent: Determines if trigger condition is persistent or not
    :type persistent: bool
    """
    def __init__(self, expression=None, initial_value = False, persistent = False):

        if isinstance(expression, str):
            self.expression = expression
        else:
            raise EventError('EventTrigger expression must be a string')

        if isinstance(initial_value, bool):
            self.value = initial_value
        else:
            raise EventError('EventTrigger initial_value must be bool')

        if isinstance(persistent, bool):
            self.persistent = persistent
        else:
            raise EventError('EventTrigger.persistent must be bool')
    def __str__(self):
        return self.expression

    def sanitized_expression(self, species_mappings, parameter_mappings):
        '''
        Sanitize the event trigger expression.

        :param species_mappings: Mapping of species names to sanitized species names.
        :type species_mappings: dict

        :param parameter_mappings: Mapping of parameter names to sanitized parameter names.
        :type parameter_mappings: dict

        :returns: The sanitized expression.
        :rtype: str
        '''
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_expression = self.expression
        for i, name in enumerate(names):
            sanitized_expression = sanitized_expression.replace(name, "{"+str(i)+"}")
        return sanitized_expression.format(*replacements)

class Event(Jsonify):
    """
    An Event can be given as an assignment_expression (function) or directly
    as a value (scalar). If given an assignment_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    :param name: The name by which this Event is called or referenced in reactions.
    :type name: str

    :param assignments: List of EventAssignments to be executed at trigger or delay
    :type assignments: str

    :param trigger: contains math expression which can be evaluated to
        a boolean result.  Upon the transition from 'False' to 'True',
        event assignments may be executed immediately, or after a
        designated delay.
    :type trigger: EventTrigger

    :param delay: contains math expression evaluable within model namespace.
        This expression designates a delay between the trigger of
        an event and the execution of its assignments.
    :type delay: str

    :param priority: Contains math expression evaluable within model namespace.
    :type priority: str

    :type use_values_from_trigger_time: bool
    """

    def __eq__(self, other):
        return str(self) == str(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __lt__(self, other):
        return str(self) < str(other)

    def __le__(self, other):
        return str(self) <= str(other)

    def __hash__(self):
        if hasattr(self, '_hash'):
            return self._hash
        if hasattr(self, 'id'):
            self._hash = hash(self.id)
        elif hasattr(self, 'name'):
            self._hash = hash(self.name)
        else:
            self._hash = hash(self)
        return self._hash

    def __init__(self, name="", delay=None, assignments=[], priority="0", trigger=None,
                 use_values_from_trigger_time=False):

        # Events can contain any number of assignments
        self.assignments = []

        # Name
        if isinstance(name, str):
            self.name = name
        else:
            raise EventError(
             'name must be a valid string')

        # Trigger
        if hasattr(trigger, 'expression'):
            self.trigger = trigger
        else:
            raise EventError(
             'trigger must be set to a valid EventTrigger')

        # Delay
        if delay is None or isinstance(delay, str):
            self.delay = delay
        else:
            raise EventError(
             'delay must be a valid string or None')

        # Priority
        self.priority = priority

        # Assignments
        if isinstance(assignments, list):
            for assign in assignments:
                if hasattr(assign, 'variable'):
                    self.assignments.append(assign)
                else:
                    raise EventError('assignment list contains an item is not an EventAssignment.')
        elif hasattr(assignments, 'variable'):
            self.assignments.append(assignments)
        else:
            raise EventError(
                'assignments must contain only EventAssignments '
                'or a list of EventAssignments')
        # Use Values from Trigger Time
        if isinstance(use_values_from_trigger_time, bool):
            self.use_values_from_trigger_time = use_values_from_trigger_time
        else:
            raise EventError(
                'use_values_from_trigger_time requires bool')

    def __str__(self):
        print_string = self.name
        print_string += '\n\tTrigger: ' + str(self.trigger)
        if len(self.assignments) > 0:
            print_string += '\n\tAssignments:'
            for assign in self.assignments:
                if isinstance(assign.variable, str):
                    print_string += '\n\t\t' + assign.variable + ': ' + assign.expression
                else:
                    print_string += '\n\t\t' + assign.variable.name + ': ' + assign.expression
        return print_string

    def add_assignment(self, assignment):
        """
        Adds an EventAssignment or a list of EventAssignment.

        :param assignment: The event or list of events to be added to this event.
        :type assignment: EventAssignment or list[EventAssignment]
        """
        if isinstance(assignment, list):
            for assign in assignment:
                self.add_assignment(assign)
        elif hasattr(assignment, 'variable'):
            self.assignments.append(assignment)
        else:
            raise EventError(
                "Unexpected parameter for add_assignment. Assignment must be an EventAssignment or list of "
                "EventAssignments objects"
            )
        return assignment
