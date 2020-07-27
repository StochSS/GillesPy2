from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.gillespyError import *
import numpy as np
import uuid
import ast


class Reaction(SortableObject):
    """
    Models a single reaction. A reaction has its own dicts of species
    (reactants and products) and parameters. The reaction's propensity
    function needs to be evaluable (and result in a non-negative scalar
    value) in the namespace defined by the union of those dicts.

    :param name: The name by which the reaction is called (optional).
    :type name: str
    :param reactants: The reactants that are consumed in the reaction, with stoichiometry. An
    example would be {R1 : 1, R2 : 2} if the reaction consumes two of R1 and
    one of R2, where R1 and R2 are Species objects.
    :type reactants: dict
    :param products: The species that are created by the reaction event, with stoichiometry. Same format as reactants.
    :type products: dict
    :param propensity_function: The custom propensity function for the reaction. Must be evaluable in the
    namespace of the reaction using C operations.
    :type propensity_function: str
    :param massaction: The switch to use a mass-action reaction. If set to True, a rate value is required.
    :type massaction: bool
    :param rate: The rate of the mass-action reaction, take care to note the units.
    :type rate: float
    :param annotation: An optional note about the reaction.
    :type annotation: str

    Notes
    ----------
    For a species that is NOT consumed in the reaction but is part of a mass
    action reaction, add it as both a reactant and a product.

    Mass-action reactions must also have a rate term added. Note that the input
    rate represents the mass-action constant rate independent of volume.

    if a name is not provided for reactions, the name will be populated by the
    model based on the order it was added. This could impact seeded simulation 
    results if the order of addition is not preserved.
    """

    def __init__(self, name="", reactants={}, products={}, propensity_function=None, massaction=False, rate=None,
                 annotation=None):
        """
        Initializes the reaction using short-hand notation.
        """

        # Metadata
        self.name = name
        self.annotation = ""

        # We might use this flag in the future to automatically generate
        # the propensity function if set to True.
        if propensity_function is not None:
            self.massaction = False
            self.marate = None
        else:
            self.massaction = True

        self.propensity_function = propensity_function
        self.ode_propensity_function = propensity_function

        if self.propensity_function is not None and self.massaction:
            errmsg = ("Reaction {} You cannot set the propensity type to mass-action and simultaneously set a "
                      "propensity function.").format(self.name)
            raise ReactionError(errmsg)

        self.reactants = {}
        for r in reactants:
            rtype = type(r).__name__
            if rtype == 'instance':
                self.reactants[r.name] = reactants[r]
            else:
                self.reactants[r] = reactants[r]

        self.products = {}
        for p in products:
            rtype = type(p).__name__
            if rtype == 'instance':
                self.products[p.name] = products[p]
            else:
                self.products[p] = products[p]

        if self.massaction:
            self.type = "mass-action"
            if rate is None:
                self.marate = None
            else:
                self.marate = rate
                self.__create_mass_action()
        else:
            self.type = "customized"

            def __customPropParser():
                pow_func = ast.parse("pow", mode="eval").body

                class ExpressionParser(ast.NodeTransformer):
                    def visit_BinOp(self, node):
                        node.left = self.visit(node.left)
                        node.right = self.visit(node.right)
                        if isinstance(node.op, (ast.BitXor, ast.Pow)):
                            # ast.Call calls defined function, args include which nodes
                            # are effected by function call
                            call = ast.Call(func=pow_func,
                                            args=[node.left, node.right],
                                            keywords=[])
                            # Copy_location copies lineno and coloffset attributes
                            # from old node to new node. ast.copy_location(new_node,old_node)
                            call = ast.copy_location(call, node)
                            # Return changed node
                            return call
                        # No modification to node, classes extending NodeTransformer methods
                        # Always return node or value
                        else:
                            return node

                    def visit_Name(self, node):
                        # Visits Name nodes, if the name nodes "id" value is 'e', replace with numerical constant
                        if node.id == 'e':
                            nameToConstant = ast.copy_location(ast.Num(float(np.e), ctx=node.ctx), node)
                            return nameToConstant
                        return node

                expr = self.propensity_function
                expr = expr.replace('^', '**')
                expr = ast.parse(expr, mode='eval')
                expr = ExpressionParser().visit(expr)

                class ToString(ast.NodeVisitor):
                    def __init__(self):
                        self.string = ''

                    def _string_changer(self, addition):
                        self.string += addition

                    def visit_BinOp(self, node):
                        self._string_changer('(')
                        self.visit(node.left)
                        self.visit(node.op)
                        self.visit(node.right)
                        self._string_changer(')')

                    def visit_Name(self, node):
                        self._string_changer(node.id)
                        self.generic_visit(node)

                    def visit_Num(self, node):
                        self._string_changer(str(node.n))
                        self.generic_visit(node)

                    def visit_Call(self, node):
                        self._string_changer(node.func.id + '(')
                        counter = 0
                        for arg in node.args:
                            self.visit(arg)
                            if counter == 0:
                                self._string_changer(',')
                                counter += 1
                        self._string_changer(')')

                    def visit_Add(self, node):
                        self._string_changer('+')
                        self.generic_visit(node)

                    def visit_Div(self, node):
                        self._string_changer('/')
                        self.generic_visit(node)

                    def visit_Mult(self, node):
                        self._string_changer('*')
                        self.generic_visit(node)

                    def visit_UnaryOp(self, node):
                        self._string_changer('(')
                        self.visit_Usub(node)
                        self._string_changer(')')

                    def visit_Sub(self, node):
                        self._string_changer('-')
                        self.generic_visit(node)

                    def visit_Usub(self, node):
                        self._string_changer('-')
                        self.generic_visit(node)

                newFunc = ToString()
                newFunc.visit(expr)
                return newFunc.string
            self.propensity_function = __customPropParser()

    def __str__(self):
        print_string = self.name
        if len(self.reactants):
            print_string += '\n\tReactants'
            for r, stoich in self.reactants.items():
                try:
                    print_string += '\n\t\t' + r.name + ': ' + str(stoich)
                except Exception as e:
                    print_string += '\n\t\t' + r + ': ' + 'INVALID - ' + str(e)
        if len(self.products):
            print_string += '\n\tProducts'
            for p, stoich in self.products.items():
                try:
                    print_string += '\n\t\t' + p.name + ': ' + str(stoich)
                except Exception as e:
                    print_string += '\n\t\t' + p + ': ' + 'INVALID - ' + str(e)
        print_string += '\n\tPropensity Function: ' + self.propensity_function
        return print_string

    def verify(self):
        """ Check if the reaction is properly formatted.
        Does nothing on sucesss, raises and error on failure."""
        if self.marate is None and self.propensity_function is None:
            raise ReactionError("You must specify either a mass-action rate or a propensity function")
        if len(self.reactants) == 0 and len(self.products) == 0:
            raise ReactionError("You must have a non-zero number of reactants or products.")

    def __create_mass_action(self):
        """
        Initializes the mass action propensity function given
        self.reactants and a single parameter value.
        """
        # We support zeroth, first and second order propensities only.
        # There is no theoretical justification for higher order propensities.
        # Users can still create such propensities if they really want to,
        # but should then use a custom propensity.
        total_stoch = 0
        for r in sorted(self.reactants):
            total_stoch += self.reactants[r]
        if total_stoch > 2:
            raise ReactionError("Reaction: A mass-action reaction cannot involve more than two of one species or one "
                                "of two species. To declare a custom propensity, replace 'rate' with "
                                "'propensity_function'.")

        # Case EmptySet -> Y

        propensity_function = self.marate.name
        ode_propensity_function = self.marate.name

        # There are only three ways to get 'total_stoch==2':
        for r in sorted(self.reactants):
            if isinstance(r, str):
                rname = r
            else:
                rname = r.name
            # Case 1: 2X -> Y
            if self.reactants[r] == 2:
                propensity_function = (propensity_function +
                                       "*" + rname + "*(" + rname + "-1)/vol")
                ode_propensity_function += '*' + rname + '*' + rname
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += "*" + rname
                ode_propensity_function += '*' + rname

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += "/vol"
        elif order == 0:
            propensity_function += "*vol"

        self.propensity_function = propensity_function
        self.ode_propensity_function = ode_propensity_function

    def setType(self, rxntype):
        """
        Sets reaction type to either "mass-action" or "customized"

        :param rxntype: Either "mass-action" or "customized"
        :type rxntype: str
        """
        if rxntype.lower() not in {'mass-action', 'customized'}:
            raise ReactionError("Invalid reaction type.")
        self.type = rxntype.lower()

        self.massaction = False if self.type == 'customized' else True

    def addReactant(self, S, stoichiometry):
        """
        Adds a reactant to the reaction (species that is consumed)

        :param S: Reactant to add to this reaction.
        :type S: gillespy2.Species
        :param stoichiometry: The stoichiometry of the given reactant.
        :type stoichiometry: int
        """
        if stoichiometry <= 0:
            raise ReactionError("Reaction Stoichiometry must be a \
                                    positive integer.")
        self.reactants[S.name] = stoichiometry

    def addProduct(self, S, stoichiometry):
        """
        Adds a product to the reaction (species that is created)

        :param S: Product to add to this reaction.
        :type S: gillespy2.Species
        :param stoichiometry: The stoichiometry of the given product.
        :type stoichiometry: int
        """
        self.products[S.name] = stoichiometry

    def Annotate(self, annotation):
        """
        Adds a note to the reaction

        :param annotation: An optional note about the reaction.
        :type annotation: str
        """
        self.annotation = annotation

    def sanitized_propensity_function(self, species_mappings, parameter_mappings):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_propensity = self.propensity_function
        for id, name in enumerate(names):
            sanitized_propensity = sanitized_propensity.replace(name, "{" + str(id) + "}")
        return sanitized_propensity.format(*replacements)
