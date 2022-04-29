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

from typing import (
    List,
    Union
)

from gillespy2.core.species import Species
from gillespy2.core.parameter import Parameter

from gillespy2.core.gillespyError import *
from gillespy2.core.jsonify import Jsonify


class EventAssignment(Jsonify):
    """
    An :class:`EventAssignment` describes a change to the current :class:`~gillespy2.core.model.Model` that will
    be conditionally executed during a simulation. This assignment can either be fired at the time its associated
    :class:`EventTrigger` transitions from :code:`False` to :code:`True`, or after a delay if one is configured within 
    the :class:`Event` this object is assigned to.

    :param name: The name by which this :class:`EventAssignment` object will be identified.

    :param variable: The target model component that will be modified by :code:`expression`.

    :param expression: The string to be evaluated when the event is fired. This expression must be evaluatable
        within the context of a :class:`~gillespy2.core.model.Model` namespace. Results of this expression's
        evaluation will be assigned to the :class:`EventAssignment` variable.
    """

    def __init__(
        self, 
        name: str = None, 
        variable: Union[Species, Parameter] = None, 
        expression: str = None
    ):

        if name in (None, ""):
            self.name = f'evn{uuid.uuid4()}'.replace('-', '_')
        else:
            self.name = name

        self.variable = variable
        self.expression = expression

        if expression is not None:
            self.expression = str(expression)

        from gillespy2.core.parameter import Parameter
        from gillespy2.core.species import Species
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

    def __str__(self) ->  str:
        return f"{self.variable}: {self.expression}"

class EventTrigger(Jsonify):
    """
    An :class:`EventTrigger` reacts to changes in the :class:`gillespy2.core.model.Model` simulation environment to 
    fire one or more :class:`EventAssignment` objects.

    Trigger detects changes in model/environment conditions in order to fire an
    event.  A Trigger contains an expression, a mathematical function which can
    be evaluated to a boolean value within a model's namespace.  Upon
    transitioning from 'false' to 'true', this trigger will cause the immediate
    execution of an event's list of assignments if no delay is present, otherwise,
    the delay evaluation will be initialized.

    :param expression: A math expression which evaluates to a boolean result with respect to the simulation context
        of the running model. Upon transition from :code:`False` to :code:`True`, one or more :class:`EventAssignment`
        objects can either be executed immediately or after a designated delay. This expression should be evaluatable
        within the namespace of the current :class:`~gillespy2.core.model.Model`. 

    :param value: The initial boolean state of this :class:`EventTrigger` at time :code:`t = 0`. 

    :param persistent: Determines if the trigger condition is persistent.
    """
    def __init__(
        self, 
        expression: str = None, 
        initial_value: bool = False, 
        persistent: bool = False
    ):

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
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_expression = self.expression
        for id, name in enumerate(names):
            sanitized_expression = sanitized_expression.replace(name, "{"+str(id)+"}")
        return sanitized_expression.format(*replacements)


class Event(Jsonify):
    """
    An Event can be given as an assignment_expression (function) or directly
    as a value (scalar). If given an assignment_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    :param name: The name by which this :class:`Event` is called or referenced in reactions.

    :param assignments: A list of :class:`EventAssignment` objects which, upon execution, mutate the simulation state
        of the current :class:`~gillespy2.core.model.Model`. Assignments can be executed immediately or after a delay
        as defined by the :code:`delay` argument.

    :param trigger: A :class:`EventTrigger` which evaluates, in the context of the current 
        :class:`~gillespy2.core.model.Model` namespace, to a boolean result. Upon transition from :code:`False` to 
        :code:`True`, one or more :class:`EventAssignment` objects will be executed immediately or after a specified 
        delay.

    :param delay: A math expression which evaluates the delay between the state transition of an :class:`EventTrigger`
        and the subsequent execution of one or more :class:`EventAssignment` objects. Must be evaluatable within the
        context of the current :class:`~gillespy2.core.model.Model` namespace.

    :param priority: A math expression that is evaluatable within the current :class:`~gillespy2.core.model.Model`
        namespace.
    """

    def __init__(
        self, 
        name: str = "", 
        delay: str = None, 
        assignments: List[EventAssignment] = [], 
        priority: str = "0", 
        trigger: EventTrigger = None,
        use_values_from_trigger_time: bool = False
    ):
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


    def __eq__(self, other: EventAssignment) -> bool:
        return str(self) == str(other)

    def __ne__(self, other: EventAssignment) -> bool:
        return not self.__eq__(other)

    def __gt__(self, other: EventAssignment) -> bool:
        return not self.__le__(other)

    def __ge__(self, other: EventAssignment) -> bool:
        return not self.__lt__(other)

    def __lt__(self, other: EventAssignment) -> bool:
        return str(self) < str(other)

    def __le__(self, other: EventAssignment) -> bool:
        return str(self) <= str(other)

    def __hash__(self) -> int:
        if hasattr(self, '_hash'):
            return self._hash
        if hasattr(self, 'id'):
            self._hash = hash(self.id)
        elif hasattr(self, 'name'):
            self._hash = hash(self.name)
        else:
            self._hash = hash(self)
        return self._hash

    def __str__(self) -> str:
        print_string = self.name
        print_string += '\n\tTrigger: ' + str(self.trigger)
        if len(self.assignments):
            print_string += '\n\tAssignments:'
            for a in self.assignments:
                if isinstance(a.variable, str):
                    print_string += '\n\t\t' + a.variable + ': ' + a.expression
                else:
                    print_string += '\n\t\t' + a.variable.name + ': ' + a.expression
        return print_string

    def add_assignment(
        self, 
        assignment: Union[EventAssignment, List[EventAssignment]]
    ) -> Union[EventAssignment, List[EventAssignment]]:
        """
        Adds one or more :class:`EventAssignment` objects.

        :param assignment: The :class:`EventAssignment` or list of :class:`EventAssignment` objects to be added to
            this event.
        """

        if hasattr(assignment, 'variable'):
            self.assignments.append(assignment)
        elif isinstance(assignment, list):
            for assign in assignment:
                if hasattr(assign, 'variable'):
                    self.assignments.append(assign)
                else:
                    raise EventError('add_assignment failed to add EventAssignment. Assignment to be added must be of '
                                     'type EventAssignment or list of EventAssignment objects.')
        else:
            raise ModelError("Unexpected parameter for add_assignment. Parameter must be EventAssignment or list of "
                             "EventAssignments")
        return assignment
