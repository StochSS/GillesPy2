from gillespy2.core.sortableobject import SortableObject

class FunctionDefinition(SortableObject):
    """
    Object representation defining an evaluable function to be used during
    simulation of a GillesPy2 model

    :param name: Name of the function to be made and called
    :type name: str
    :param function: Defined function body of operation to be performed.
    :type function: str
    :param variables: String names of Variables to be used as arguments to function.
    :type variables: list of strings
    """

    def __init__(self, name="", function=None, args=[]):

        import math
        eval_globals = math.__dict__

        self.name = name
        args = ', '.join(args)
        self.function = eval('lambda ' + args + ': ' + function, eval_globals)
        if self.function is None:
            raise TypeError

    def sanitized_function(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_function = self.function
        for id, name in enumerate(names):
            sanitized_function = sanitized_function.replace(name, "{" + str(id) + "}")
        return sanitized_function.format(*replacements)