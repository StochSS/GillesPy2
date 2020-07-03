from gillespy2.core.sortableobject import SortableObject


class AssignmentRule(SortableObject):
    """
    An AssignmentRule is used to express equations that set the values of
    variables.  This would correspond to a function in the form of x = f(V)
    :param name: Name of the Rule
    :type name: str
    :param variable: Target Species/Parameter to be modified by rule
    :type variable: str
    :param formula: String representation of formula to be evaluated
    :type formula: str
    """

    def __init__(self, variable=None, formula=None, name=None):
        self.variable = variable
        self.formula = formula
        self.name = name

    def __str__(self):
        return self.variable + ': ' + self.formula

    def sanitized_formula(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for id, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(id) + "}")
        return sanitized_formula.format(*replacements)