from collections import OrderedDict
from gillespy2.core.gillespyError import *
from gillespy2.core.gillespy2 import *

class EventAssignment:
    """
    An EventAssignment describes a change to be performed to the current model
    simulation.  This is assignment can either be fired at the time its
    associated trigger changes from false to true, or after a specified delay,
    depending on how the Event to which it is assigned is configured.

    Attributes
    ----------
    variable : gillespy2.Species, gillespy2.Parameter
        Target model component to be modified by the EventAssignment
        expression. Valid target variables include gillespy2 Species,
        Parameters, and Compartments.
    expression : str
        String to be evaluated when the event is fired.  This expression must
        be evaluable within the model namespace, and the results of it's
        evaluation will be assigned to the EventAssignment variable.
    """

    def __init__(self, variable=None, expression=None):

        self.variable = variable
        self.expression = expression

        if expression is not None:
            self.expression = str(expression)

        
        from gillespy2.core.gillespy2 import Species, Parameter
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

    Attributes
    ----------
    expression : str
        String for a function calculating EventTrigger values. Should be evaluable
        in namespace of Model.
    value : bool
        Value of EventTrigger at simulation start, with time t=0
    persistent : bool
        Determines of trigger condition is persistent or not.
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
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key = lambda x: len(x), reverse=True)
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

    Attributes
    ----------
    name : str
        The name by which this Event is called or referenced in reactions.
    assignments : list
        List of EventAssignments to be executed at trigger or delay
    trigger : EventTrigger
        contains math expression which can be evaluated to 
        a boolean result.  Upon the transition from 'False' to 'True',
        event assignments may be executed immediately, or after a
        designated delay.
    delay : string
        contains math expression evaluable within model namespace.
        This expression designates a delay between the trigger of
        an event and the execution of its assignments.
    priority : string
        contains math expression evaluable within model namespace.
        TODO: MORE INFO
    use_values_from_trigger_time: boolean
    """

    def __init__(self, name="", delay = None, assignments = [], priority="0", 
        trigger = None, use_values_from_trigger_time = False):

        # Events can contain any number of assignments
        self.assignments = []

        # Name
        if isinstance(name, str):
            self.name = name
        else:
            raise EventError(
             'name must be a valid string')
        
        # Trigger
        if isinstance(trigger, EventTrigger):
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
                if isinstance(assign, EventAssignment):
                    self.assignments.append(assign)
                else:
                    raise EventError('assignment list contains an item '
                        'is not an EventAssignment.')
        elif isinstance(assignments, EventAssignment):
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
                print_string += '\n\t\t' + a.variable.name + ': ' + a.expression
        return print_string

    def add_assignment(self, assignment):
        """
        Adds an eventAssignment or a list of eventAssignments.

        Attributes
        ----------
        assignment : EventAssignment or a list of EventAssignments
            The event or list of events to be added to this event.
        """

        if isinstance(assignment, EventAssignment):
            self.assignments.append(assignment)
        elif isinstance(assignment, list):
            for assign in assignment:
                if isinstance(assign, EventAssignment):
                    self.assignments.append(assign)
                else:
                    raise EventError('add_assignment failed to add EventAssignment. '
                    'Assignment to be added must be of type EventAssignment '
                    'or list of EventAssignment objects.')
        else:
            raise ModelError("Unexpected parameter for add_assignment. Parameter must be EventAssignment or list of EventAssignments")
        return obj



