from gillespy2.core.sortableobject import SortableObject


class RateRule(SortableObject):
    """
    A RateRule is used to express equations that determine the rates of change
    of variables. This would correspond to a function in the form of dx/dt=f(W)

    :param name: Name of Rule
    :type name: str
    :param variable: Target Species/Parameter to be modified by rule
    :type variable: str
    :param formula: String representation of formula to be evaluated
    :type formula: str
    """

    def __init__(self, variable=None, formula='', name=None):
        self.formula = formula
        self.variable = variable
        self.name = name

    def __str__(self):
        try:
            return self.name + ': Var: ' + self.variable + ': ' + self.formula
        except:
            return 'Rate Rule: {} contains an invalid variable or formula'.format(self.name)

    def sanitized_formula(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_formula = self.formula
        for id, name in enumerate(names):
            sanitized_formula = sanitized_formula.replace(name, "{" + str(id) + "}")
        return sanitized_formula.format(*replacements)