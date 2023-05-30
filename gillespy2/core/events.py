# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

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

from gillespy2.core.species import Species
from gillespy2.core.jsonify import Jsonify

from gillespy2.core.gillespyError import EventError, EventAssignmentError, EventTriggerError

class EventAssignment(Jsonify):
    """
    An EventAssignment describes a change to be performed to the current model
    simulation.  This is assignment can either be fired at the time its
    associated trigger changes from false to true, or after a specified delay,
    depending on how the Event to which it is assigned is configured.

    :param variable: Target gillespy2 Species to be modified by the EventAssignment expression.
    :type variable: gillespy2.Species, gillespy2.Parameter

    :param expression: String to be evaluated when the event is fired.  This expression must
        be evaluable within the model namespace, and the results of it's
        evaluation will be assigned to the EventAssignment variable.
    :type expression: str

    """
    def __init__(self, name=None, variable=None, expression=None):
        if name in (None, ""):
            name = f'evn_assign{uuid.uuid4()}'.replace('-', '_')
        else:
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning("EventAssignment.name has been deprecated.")
        self.__name_deprecated = name
        if isinstance(expression, (int, float)):
            expression = str(expression)

        self.expression = expression

        self.validate(variable=variable)

        if variable is not None:
            vtype = type(variable).__name__
            if vtype == 'Species':
                variable = variable.name
        self.variable = variable

    def __str__(self):
        return f"{self.variable}: {self.expression}"

    def __getattr__(self, key):
        if key == 'name':
            from gillespy2.core import log # pylint: disable=import-outside-toplevel
            log.warning('EventAssignment.name has been deprecated.')
            return self.__name_deprecated
        raise AttributeError

    def _create_sanitized_event_assignment(self, species_mappings, parameter_mappings):
        variable = species_mappings[self.variable.name]
        expression = self.sanitized_expression(species_mappings, parameter_mappings)
        return EventAssignment(expression=expression, variable=variable)

    def sanitized_expression(self, species_mappings, parameter_mappings):
        '''
        Sanitize the event assignment expression.

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

    def validate(self, variable=None, coverage="all"):
        """
        Validate the event assignment.

        :param variable: Target Species to be modified by event assignment
        :type variable: str

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises EventAssignmentError: Attribute is of invalid type.  Required attribute set to None.  \
                              Attribute is value outside of accepted bounds.
        """
        # Check variable
        if coverage in ("all", "variable"):
            if variable is None:
                if not hasattr(self, "variable") or self.variable is None:
                    raise EventAssignmentError("Event assignments must have a variable.")
                variable = self.variable

            if not (isinstance(variable, (str, Species)) or type(variable).__name__ == 'Species'):
                raise EventAssignmentError("variable must be of type str or GillesPy2.Species.")
            if variable == "":
                raise EventAssignmentError("variable can't be an empty string.")

        # Check expression
        if coverage in ("all", "expression"):
            if not isinstance(self.expression, str):
                raise EventAssignmentError("expression must be of type str.")
            if self.expression == "":
                raise EventAssignmentError("expression can't be an empty string.")

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
    def __init__(self, expression=None, initial_value=False, persistent=False):
        if isinstance(expression, (int, float)):
            expression = str(expression)

        self.expression = expression
        self.value = initial_value
        self.persistent = persistent

        self.validate()

    def __str__(self):
        return f"Value: {self.value}, Persistent: {self.persistent},\n\t\tExpression: {self.expression}"

    def _create_sanitized_event_trigger(self, species_mappings, parameter_mappings):
        expression = self.sanitized_expression(species_mappings, parameter_mappings)
        return EventTrigger(expression=expression, initial_value=self.value, persistent=self.persistent)

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

    def validate(self, coverage="all"):
        """
        Validate the event trigger.

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises EventTriggerError: Attribute is of invalid type.  Required attribute set to None.  \
                              Attribute is value outside of accepted bounds.
        """
        # Check expression
        if coverage in ("all", "expression"):
            if not isinstance(self.expression, str):
                raise EventTriggerError("expression must be of type str.")
            if self.expression == "":
                raise EventTriggerError("expression can't be an empty string.")

        # Check value
        if coverage in ("all", "value", "initial_value"):
            if not isinstance(self.value, bool):
                raise EventTriggerError(f"value must be of type bool not {type(self.value)}.")

        # Check persistent
        if coverage in ("all", "persistent"):
            if not isinstance(self.persistent, bool):
                raise EventTriggerError(f"value must be of type bool not {type(self.persistent)}.")

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
    def __init__(self, name=None, delay=None, assignments=None, priority="0", trigger=None,
                 use_values_from_trigger_time=False):
        if name in (None, ""):
            name = f'evn{uuid.uuid4()}'.replace('-', '_')
        if assignments is None:
            assignments = []
        elif not isinstance(assignments, list):
            assignments = [assignments]

        self.name = name
        self.trigger = trigger
        self.assignments = assignments
        self.delay = delay
        self.priority = priority
        self.use_values_from_trigger_time = use_values_from_trigger_time

        self.validate()

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

    def __str__(self):
        print_string = self.name
        print_string += '\n\tTrigger:\n\t\t' + str(self.trigger)
        if len(self.assignments) > 0:
            print_string += '\n\tAssignments:'
            for assign in self.assignments:
                if isinstance(assign.variable, str):
                    print_string += '\n\t\t' + assign.variable + ': ' + assign.expression
                else:
                    print_string += '\n\t\t' + assign.variable.name + ': ' + assign.expression
        return print_string

    def _create_sanitized_event(self, n_ndx, species_mappings, parameter_mappings):
        name = f"E{n_ndx}"
        priority, delay = self.sanitized_expression(species_mappings, parameter_mappings)
        trigger = self.trigger._create_sanitized_event_trigger(species_mappings, parameter_mappings)
        assignments = [
            assignment._create_sanitized_event_assignment(
                species_mappings, parameter_mappings
            ) for assignment in self.assignments
        ]
        return Event(
            name=name, priority=priority, delay=delay, trigger=trigger, assignments=assignments,
            use_values_from_trigger_time=self.use_values_from_trigger_time
        )

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

    def sanitized_expression(self, species_mappings, parameter_mappings):
        '''
        Sanitize the event's delay and priority expressions.

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
        if self.delay is None:
            sanitized_delay = None
        else:
            sanitized_delay = self.delay
            for i, name in enumerate(names):
                sanitized_delay = sanitized_delay.replace(name, "{"+str(i)+"}")
            sanitized_delay.format(*replacements)

        sanitized_priority = self.priority
        for i, name in enumerate(names):
            sanitized_priority = sanitized_priority.replace(name, "{"+str(i)+"}")
        return sanitized_priority.format(*replacements), sanitized_delay

    def validate(self, coverage="all"):
        """
        Validate the event.

        :param coverage: The scope of attributes to validate.  Set to an attribute name to restrict validation \
                         to a specific attribute.
        :type coverage: str

        :raises EventTriggerError: Attribute is of invalid type.  Required attribute set to None.  \
                              Attribute is value outside of accepted bounds.
        """
        # Check name
        if coverage in ("all", "name"):
            if self.name is None:
                raise EventError("name can't be None type.")
            if not isinstance(self.name, str):
                raise EventError(f"name must be of type str not {type(self.name)}.")
            if self.name == "":
                raise EventError("name can't be an empty string.")

        # Check trigger
        if coverage in ("all", "trigger"):
            if self.trigger is None:
                raise EventError("trigger can't be None type.")
            if not (isinstance(self.trigger, EventTrigger) or type(self.trigger).__name__ == "EventTrigger"):
                raise EventError(f"trigger must be of type gillespy2.EventTrigger not {type(self.trigger)}.")
            try:
                self.trigger.validate()
            except EventTriggerError as err:
                raise EventError(f"trigger must be a valid gillespy2.EventTrigger: {str(err)}") from err

        # Check assignments
        if coverage in ("all", "assignments"):
            for assignment in self.assignments:
                if assignment is None:
                    raise EventError("event assignments can't be None type.")
                if not (isinstance(assignment, EventAssignment) or type(assignment).__name__ == "EventAssignment"):
                    raise EventError(
                        f"event assignment must be of type gillespy2.EventAssignment not {type(assignment)}."
                    )
                try:
                    assignment.validate()
                except EventAssignmentError as err:
                    raise EventError(f"event assignment must be a valid gillespy2.EventAssignment: {str(err)}") from err

        # Check delay
        if coverage in ("all", "delay") and self.delay is not None:
            if not isinstance(self.delay, str):
                raise EventAssignmentError("delay must be of type str.")
            if self.delay == "":
                raise EventAssignmentError("delay can't be an empty string.")

        # Check priority
        if coverage in ("all", "priority"):
            if not isinstance(self.priority, str):
                raise EventAssignmentError("priority must be of type str.")
            if self.priority == "":
                raise EventAssignmentError("priority can't be an empty string.")

        # Check use_values_from_trigger_time
        if coverage in ("all", "use_values_from_trigger_time"):
            if not isinstance(self.use_values_from_trigger_time, bool):
                raise EventError(
                    f"use_values_from_trigger_time must be of type bool not {type(self.use_values_from_trigger_time)}."
                )
