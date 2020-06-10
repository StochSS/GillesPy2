from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.gillespyError import *

class Parameter(SortableObject):
    """
    A parameter can be given as an expression (function) or directly
    as a value (scalar). If given an expression, it should be
    understood as evaluable in the namespace of a parent Model.

    Attributes
    ----------
    name : str
        The name by which this parameter is called or referenced in reactions.
    expression : str
        String for a function calculating parameter values. Should be evaluable
        in namespace of Model.
    value : float
        Value of a parameter if it is not dependent on other Model entities.
    """

    def __init__(self, name="", expression=None, value=None):

        self.name = name
        # We allow expression to be passed in as a non-string type. Invalid strings
        # will be caught below. It is perfectly fine to give a scalar value as the expression.
        # This can then be evaluated in an empty namespace to the scalar value.
        self.expression = expression
        if expression is not None:
            self.expression = str(expression)

        self.value = value

        # self.value is allowed to be None, but not self.expression. self.value
        # might not be evaluable in the namespace of this parameter, but defined
        # in the context of a model or reaction.
        if self.expression is None:
            raise TypeError

        if self.value is None:
            self.evaluate()

    def __str__(self):
        return self.name + ': ' + self.expression

    def evaluate(self, namespace={}):
        """
        Evaluate the expression and return the (scalar) value in the given
        namespace.

        Attributes
        ----------
        namespace : dict (optional)
            The namespace in which to test evaluation of the parameter, if it
            involves other parameters, etc.
        """
        try:
            self.value = (float(eval(self.expression, namespace)))
        except:
            self.value = None

    def set_expression(self, expression):
        """
        Sets the expression for a parameter.
        """
        self.expression = expression
        # We allow expression to be passed in as a non-string type. Invalid
        # strings will be caught below. It is perfectly fine to give a scalar
        # value as the expression. This can then be evaluated in an empty
        # namespace to the scalar value.
        if expression is not None:
            self.expression = str(expression)

        if self.expression is None:
            raise TypeError

        self.evaluate()