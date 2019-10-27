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

        #TODO: ADD Compartment to valid variable types once implemented
        valid_variable_types = [Species, Parameter]

        if not type(self.variable) in valid_variable_types:
            raise EventError(
                'GillesPy2 Event Assignment variable must be a '
                 'valid gillespy2 species')
        if not isinstance(self.expression, str):
            raise EventError(
                             'GillesPy2 Event Assignment expression requires a '
                             'valid string expression')


class EventDelay:
    """
    An EventDelay can be given as an assignment_expression (function) or directly
    as a value (scalar). If given an assignment_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this EventDelay is called or referenced in reactions.
    delay_expression : str
        String for a function calculating EventDelay values. Should be evaluable
        in namespace of Model.
    value : float
        Value of an EventDelay if it is not dependent on other Model entities.
    """

    def __init__(self, name="", delay_expression=None, value=None, useValuesFromTriggerTime = False):

        self.name = name
        # We allow assignment_expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the assignment_expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.delay_expression = delay_expression
        if delay_expression is not None:
            self.delay_expression = str(delay_expression)

        self.useValuesFromTriggerTime = useValuesFromTriggerTime
        self.value = value

        # self.value is allowed to be None, but not self.assignment_expression. self.value
        # might not be evaluable in the namespace of this eventDelay, but defined
        # in the context of a model or reaction.
        if self.delay_expression is None:
            raise TypeError

        if self.value is None:
            self.evaluate()



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
    initial_value : bool
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
            self.initial_value = initial_value
        else:
            raise EventError('EventTrigger initial_value must be bool')


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
        trigger = None, use_values_from_trigger_time = True):

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
        if delay is None or isinstance(delay, EventDelay):
            self.delay = delay
        else:
            raise EventError(
             'delay must be a valid EventDelay or None')

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


    def add_assignment(self, assignment):
        """
        Adds an eventAssignment or a list of eventAssignments.

        Attributes
        ----------
        assignment : EventAssignment or a list of EventAssignments
            The event or list of events to be added to this event.
        """

        if isinstance(assignment, EventAssignment):
            self.listOfAssignments.append(assignment)
        elif isinstance(assignment, list):
            for assign in assignment:
                if isinstance(assign, EventAssignment):
                    self.listOfAssignments.append(assign)
                else:
                    raise EventError('add_assignment failed to add EventAssignment. '
                    'Assignment to be added must be of type EventAssignment '
                    'or list of EventAssignment objects.')
        else:
            raise ModelError("Unexpected parameter for add_assignment. Parameter must be EventAssignment or list of EventAssignments")
        return obj



