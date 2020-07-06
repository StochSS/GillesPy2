from gillespy2.core.gillespyError import *


class EventAssignment:
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

    def __init__(self, variable=None, expression=None):

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

    def __str__(self):
        return self.variable.name + ': ' + self.expression


class EventTrigger:
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
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_expression = self.expression
        for id, name in enumerate(names):
            sanitized_expression = sanitized_expression.replace(name, "{"+str(id)+"}")
        return sanitized_expression.format(*replacements)


class Event:

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

    :param use_values_from_trigger_time
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
        if len(self.assignments):
            print_string += '\n\tAssignments:'
            for a in self.assignments:
                if isinstance(a.variable, str):
                    print_string += '\n\t\t' + a.variable + ': ' + a.expression
                else:
                    print_string += '\n\t\t' + a.variable.name + ': ' + a.expression
        return print_string

    def add_assignment(self, assignment):
        """
        Adds an eventAssignment or a list of eventAssignments.
        :param assignment: EventAssignment or a list of EventAssignments
        The event or list of events to be added to this event.
        :type assignment: EventAssignment or list of EventAssignments
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



