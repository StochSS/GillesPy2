from collections import OrderedDict
from gillespy2.core.gillespyError import *

class EventAssignment:
    """
    An EventAssignment can be given as an delay_expression (function) or directly
    as a value (scalar). If given an delay_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this EventAssignment is called or referenced in reactions.
    assignment_expression : str
        String for a function calculating EventAssignment values. Should be evaluable
        in namespace of Model.
    value : float
        Value of an EventAssignment if it is not dependent on other Model entities.
    """

    def __init__(self, name="", assignment_expression=None, value=None):

        self.name = name
        # We allow assignment_expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the assignment_expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.assignment_expression = assignment_expression
        if assignment_expression is not None:
            self.assignment_expression = str(assignment_expression)

        self.value = value

        # self.value is allowed to be None, but not self.assignment_expression. self.value
        # might not be evaluable in the namespace of this event, but defined
        # in the context of a model or reaction.
        if self.assignment_expression is None:
            raise TypeError

        if self.value is None:
            self.evaluate()

    def evaluate(self, namespace={}):
        """
        Evaluate the assignment_expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the event, if it
            involves other event, etc.
        """
        try:
            self.value = (float(eval(self.assignment_expression, namespace)))
        except:
            self.value = None

    def set_expression(self, assignment_expression):
        """
        Sets the assignment_expression for a eventAssignment.
        """
        self.assignment_expression = assignment_expression
        # We allow assignment_expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the assignment_expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if assignment_expression is not None:
            self.assignment_expression = str(assignment_expression)

        if self.assignment_expression is None:
            raise TypeError

        self.evaluate()

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

    def __init__(self, name="", delay_expression=None, value=None,useValuesFromTriggerTime = False):

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

    def evaluate(self, namespace={}):
        """
        Evaluate the assignment_expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the eventDelay, if it
            involves other eventDelay, etc.
        """
        try:
            self.value = (float(eval(self.delay_expression, namespace)))
        except:
            self.value = None

    def set_expression(self, delay_expression):
        """
        Sets the assignment_expression for a eventDelay.
        """
        self.delay_expression = delay_expression
        # We allow assignment_expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the assignment_expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if delay_expression is not None:
            self.delay_expression = str(delay_expression)

        if self.delay_expression is None:
            raise TypeError

        self.evaluate()


class EventTrigger:
    """
    An EventTrigger can be given as an assignment_expression (function) or directly
    as a value (scalar). If given an assignment_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this EventTrigger is called or referenced in reactions.
    assignment_expression : str
        String for a function calculating EventTrigger values. Should be evaluable
        in namespace of Model.
    value : float
        Value of an EventTrigger if it is not dependent on other Model entities.
    """

    def __init__(self, name="", trigger_expression=None, value=None,initial_value = False,persistent = False):

        self.name = name
        # We allow assignment_expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the assignment_expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.trigger_expression = trigger_expression
        if trigger_expression is not None:
            self.trigger_expression = str(trigger_expression)

        self.value = value

        self.initial_value = initial_value

        # self.value is allowed to be None, but not self.assignment_expression. self.value
        # might not be evaluable in the namespace of this event, but defined
        # in the context of a model or reaction.
        if self.trigger_expression is None:
            raise TypeError

        if self.value is None:
            self.evaluate()

    def evaluate(self, namespace={}):
        """
        Evaluate the assignment_expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the event, if it
            involves other event, etc.
        """
        try:
            self.value = (float(eval(self.trigger_expression, namespace)))
        except:
            self.value = None

    def set_expression(self, trigger_expression):
        """
        Sets the assignment_expression for a event.
        """
        self.trigger_expression = trigger_expression
        # We allow assignment_expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the assignment_expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if trigger_expression is not None:
            self.trigger_expression = str(trigger_expression)

        if self.trigger_expression is None:
            raise TypeError

        self.evaluate()

class Event:
    """
    An Event can be given as an assignment_expression (function) or directly
    as a value (scalar). If given an assignment_expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this Event is called or referenced in reactions.
    assignment_expression : str
        String for a function calculating Event values. Should be evaluable
        in namespace of Model.
    value : float
        Value of an Event if it is not dependent on other Model entities.
    """

    def __init__(self, name="", delay = None, listOfevent_assignments = None, priority_expression=None, value=None, event_trigger = None):

        self.name = name
        # We allow assignment_expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the assignment_expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.priority_expression = priority_expression
        if priority_expression is not None:
            self.priority_expression = str(priority_expression)
        else:
            self.priority_expression = "0"
        # self.value is allowed to be None, but not self.assignment_expression. self.value
        # might not be evaluable in the namespace of this event, but defined
        # in the context of a model or reaction.

        if event_trigger is None:
            raise TypeError

        self.value = value
        self.delay = delay
        self.event_trigger = event_trigger


        self.event_assignments = OrderedDict()

        if listOfevent_assignments is not None:
            self.add_EventAssignment(listOfevent_assignments)

        if self.value is None:
            self.evaluate()

    def evaluate(self, namespace={}):
        """
        Evaluate the assignment_expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the event, if it
            involves other event, etc.
        """
        try:
            self.value = (float(eval(self.priority_expression, namespace)))
        except:
            self.value = None

    def set_expression(self, priority_expression):
        """
        Sets the assignment_expression for a event.
        """
        self.priority_expression = priority_expression
        # We allow assignment_expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the assignment_expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if priority_expression is not None:
            self.priority_expression = str(priority_expression)

        if self.priority_expression is None:
            raise TypeError

        self.evaluate()

    def set_delay(self,delay):
        self.delay = delay

        if self.delay is None:
            raise TypeError

        self.delay

    def add_EventAssignment(self, obj):
        """
        Adds an eventAssignment or a list of eventAssignments.

        Attributes
        ----------
        obj : EventAssignment or a list of EventAssignments
            The event or list of events to be added to the model object.
        """

        if isinstance(obj, EventAssignment):
            self.event_assignments[obj.name] = obj
        elif isinstance(obj,list):
            for EA in obj:
                self.add_EventAssignment(EA)
        else:
            raise ModelError("Unexpected parameter for add_EventAssignment. Parameter must be EventAssignment or list of EventAssignments")
        return obj

    def delete_eventAssignment(self, obj):
        """
           Removes an eventAssignment object by name.

           Attributes
           ----------
           obj : str
               Name of the eventAssignment object to be removed.
           """
        self.event_assignments.pop(obj)

    def delete_all_eventAssignments(self):
        """
            Removes all events from the model object.
        """
        self.event_assignments.clear()


