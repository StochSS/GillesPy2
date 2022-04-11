"""
GillesPy2 is a modeling toolkit for biochemical simulation.
Copyright (C) 2019-2021 GillesPy2 developers.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import ast
import uuid
from json.encoder import JSONEncoder

import numpy as np

from gillespy2.core.species import Species
from gillespy2.core.jsonify import Jsonify
from gillespy2.core.sortableobject import SortableObject
from gillespy2.core.gillespyError import ReactionError


class Reaction(SortableObject, Jsonify):
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

    :param massaction: (deprecated) The switch to use a mass-action reaction. If set to True, a rate value is required.
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
    def __init__(self, name=None, reactants=None, products=None, propensity_function=None,
                 rate=None, annotation=None, massaction=None):
        if massaction is not None:
            from gillespy2.core import log
            log.warning(
                """
                massaction has been deprecated.  Future releases of GillesPy2 may not support this feature.
                Set propensity_function for massaction=False or rate for massaction=True instead.
                """
            )

        if name is None:
            name = f'rxn{uuid.uuid4()}'.replace('-', '_')
        if reactants is None:
            reactants = {}
        if products is None:
            products = {}
        if isinstance(propensity_function, (int, float)):
            propensity_function = str(propensity_function)
        if isinstance(rate, (int, float)):
            rate = str(rate)
        
        self.name = name
        
        self.reactants = {}
        for r in reactants:
            rtype = type(r).__name__
            if rtype == 'Species':
                self.reactants[r.name] = reactants[r]
            else:
                self.reactants[r] = reactants[r]
        
        self.products = {}
        for p in products:
            rtype = type(p).__name__
            if rtype == 'Species':
                self.products[p.name] = products[p]
            else:
                self.products[p] = products[p]
        
        self.propensity_function = propensity_function
        self.marate = rate
        self.ode_propensity_function = None
        self.annotation = annotation
        
        self.validate()
        
        if self.propensity_function is None:
            self.massaction = True
            self.type = "mass-action"
            self.__create_mass_action()
        else:
            self.massaction = False
            self.type = "customized"
            propensity = self.__create_custom_propensity()
            self.propensity_function = propensity
            self.ode_propensity_function = propensity

    def __str__(self):
        print_string = self.name
        if len(self.reactants):
            print_string += '\n\tReactants'
            for r, stoich in self.reactants.items():
                try:
                    if isinstance(r, str):
                        print_string += '\n\t\t' + r + ': ' + str(stoich)
                    else:
                        print_string += '\n\t\t' + r.name + ': ' + str(stoich)
                except Exception as e:
                    print_string += '\n\t\t' + r + ': ' + 'INVALID - ' + str(e)
        if len(self.products):
            print_string += '\n\tProducts'
            for p, stoich in self.products.items():
                try:
                    if isinstance(p, str):
                        print_string += '\n\t\t' + p + ': ' + str(stoich)
                    else:
                        print_string += '\n\t\t' + p.name + ': ' + str(stoich)
                except Exception as e:
                    print_string += '\n\t\t' + p + ': ' + 'INVALID - ' + str(e)
        print_string += '\n\tPropensity Function: ' + self.propensity_function
        return print_string
    
    class __ExpressionParser(ast.NodeTransformer):
        def visit_BinOp(self, node):
            node.left = self.visit(node.left)
            node.right = self.visit(node.right)
            if isinstance(node.op, (ast.BitXor, ast.Pow)):
                # ast.Call calls defined function, args include which nodes
                # are effected by function call
                pow_func = ast.parse("pow", mode="eval").body
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
    
    class __ToString(ast.NodeVisitor):
        substitutions = {
            "ln": "log",
        }

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
            func_name = ToString.substitutions.get(node.func.id) \
                if node.func.id in ToString.substitutions \
                else node.func.id
            self._string_changer(func_name + '(')
            counter = 0
            for arg in node.args:
                if counter > 0:
                    self._string_changer(',')
                self.visit(arg)
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
            raise ReactionError(
                """
                Reaction: A mass-action reaction cannot involve more than two of
                one species or one of two species. To declare a custom propensity,
                replace 'rate' with 'propensity_function'.
                """
            )

        # Case EmptySet -> Y

        rtype = type(self.marate).__name__
        if rtype == 'Parameter':
            self.marate = self.marate.name
        elif rtype in ('int', 'float'):
            self.marate = str(self.rate)

        propensity_function = self.marate
        ode_propensity_function = self.marate

        # There are only three ways to get 'total_stoch==2':
        for reactant in self.reactants:
            # Case 1: 2X -> Y
            if self.reactants[reactant] == 2:
                propensity_function = f"0.5 * {propensity_function} * {reactant} * ({reactant} - 1) / vol"
                ode_propensity_function += f" * {reactant} * {reactant}"
            else:
                # Case 3: X1, X2 -> Y;
                propensity_function += f" * {reactant}"
                ode_propensity_function += f" * {reactant}"

        # Set the volume dependency based on order.
        order = len(self.reactants)
        if order == 2:
            propensity_function += " / vol"
        elif order == 0:
            propensity_function += " * vol"

        self.propensity_function = self.__create_custom_propensity(propensity_function=propensity_function)
        self.ode_propensity_function = self.__create_custom_propensity(propensity_function=ode_propensity_function)
    
    def __create_custom_propensity(self, propensity_function=None):
        if propensity_function is None:
            propensity_function = self.propensity_function
        
        expr = propensity_function
        expr = expr.replace('^', '**')
        expr = ast.parse(expr, mode='eval')
        expr = self.__ExpressionParser().visit(expr)
        
        newFunc = self.__ToString()
        newFunc.visit(expr)
        return newFunc.string
    
    def addProduct(self, *args, **kwargs):
        """
        Add a product to this reaction (deprecated)

        :param species: Species object to be produced by the reaction
        :type species: spatialpy.core.species.Species

        :param stoichiometry: Stoichiometry of this product.
        :type stoichiometry: int
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.addProduct has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  Use Reaction.add_product instead.
            """
        )

        self.add_product(*args, **kwargs)

    def add_product(self, species, stoichiometry):
        """
        Add a product to this reaction

        :param species: Species object to be produced by the reaction
        :type species: spatialpy.core.species.Species

        :param stoichiometry: Stoichiometry of this product.
        :type stoichiometry: int
        """
        if not (isinstance(species, (str, Species)) or type(species).__name__ == 'Species'):
            raise ReactionError("species must be of type string or SpatialPy.Species. ")
        if not isinstance(stoichiometry, int) or stoichiometry <= 0:
            raise ReactionError("Stoichiometry must be a positive integer.")
        name = species if isinstance(species, str) else species.name
        self.products[name] = stoichiometry
        
    def addReactant(self, *args, **kwargs):
        """
        Add a reactant to this reaction (deprecated)

        :param species: reactant Species object
        :type species: spatialpy.core.species.Species

        :param stoichiometry: Stoichiometry of this participant reactant
        :type stoichiometry: int
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.addReactant has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  Use Reaction.add_reactant instead.
            """
        )

        self.add_reactant(*args, **kwargs)

    def add_reactant(self, species, stoichiometry):
        """
        Add a reactant to this reaction

        :param species: reactant Species object
        :type species: spatialpy.core.species.Species

        :param stoichiometry: Stoichiometry of this participant reactant
        :type stoichiometry: int
        """
        if not (isinstance(species, (str, Species)) or type(species).__name__ == 'Species'):
            raise ReactionError("Species must be of type string or spatialpy.Species. ")
        if not isinstance(stoichiometry, int) or stoichiometry <= 0:
            raise ReactionError("Stoichiometry must be a positive integer.")
        name = species if isinstance(species, str) else species.name
        self.reactants[name] = stoichiometry
    
    def Annotate(self, *args, **kwargs):
        """
        Add an annotation to this reaction (deprecated).

        :param annotation: Annotation note to be added to reaction
        :type annotation: str
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.Annotate has been deprecated.  Future releases of GillesPy2 may
            not support this feature.  Use Reaction.set_annotation instead.
            """
        )

        self.set_annotation(*args, **kwargs)

    def set_annotation(self, annotation):
        """
        Add an annotation to this reaction.

        :param annotation: Annotation note to be added to reaction
        :type annotation: str
        """
        if annotation is None:
            raise ReactionError("Annotation cannot be None.")
        self.annotation = annotation

    @classmethod
    def from_json(cls, json_object):
        new = Reaction.__new__(Reaction)
        new.__dict__ = json_object

        # Same as in to_dict(), but we need to reverse it back into its original representation.
        new.products = { x["key"]: x["value"] for x in json_object["products"] }
        new.reactants = { x["key"]: x["value"] for x in json_object["reactants"] }

        return new
    
    def sanitized_propensity_function(self, species_mappings, parameter_mappings, ode=False):
        names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
                       reverse=True)
        replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
                        for name in names]
        sanitized_propensity = self.ode_propensity_function if ode else self.propensity_function
        for id, name in enumerate(names):
            sanitized_propensity = sanitized_propensity.replace(name, "{" + str(id) + "}")
        return sanitized_propensity.format(*replacements)

    def _create_sanitized_reaction(self, n_ndx, species_mappings, parameter_mappings):
        name = f"R{n_ndx}"
        reactants = {species_mappings[species.name]: self.reactants[species] for species in self.reactants}
        products = {species_mappings[species.name]: self.products[species] for species in self.products}
        propensity_function = self.sanitized_propensity_function(species_mappings, parameter_mappings)
        return Reaction(name=name, reactants=reactants, products=products, propensity_function=propensity_function)

    def setType(self, rxntype):
        """
        Sets reaction type to either "mass-action" or "customized" (deprecated)

        :param rxntype: Either "mass-action" or "customized"
        :type rxntype: str
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.setType has been deprecated.  Future releases of GillesPy2 may not support this feature.
            Set propensity_function for Reaction.type="customized" or rate for Reaction.type="mass-action" instead.
            """
        )

        if rxntype.lower() not in {'mass-action', 'customized'}:
            raise ReactionError("Invalid reaction type.")
        self.type = rxntype.lower()

        self.massaction = False if self.type == 'customized' else True
    
    def to_dict(self):
        temp = vars(self).copy()

        # This to_dict overload is needed because, while Species is hashable (thanks to its inheretence),
        # objects are not valid key values in the JSON spec. To fix this, we set the Species equal to some key 'key', 
        # and it's value equal to some key 'value'.

        temp["products"] = list({ "key": k, "value": v} for k, v in self.products.items() )
        temp["reactants"] = list( {"key": k, "value": v} for k, v in self.reactants.items() )

        return temp
    
    def verify(self, *args, **kwargs):
        """
        Check if the reaction is properly formatted.
        Does nothing on sucesss, raises and error on failure.
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.verify has been deprecated.  Future releases of GillesPy2 may
            not support this feature. Use Reaction.validate instead.
            """
        )

        self.validate(*args, **kwargs)

    def validate(self):
        """ Check if the reaction is properly formatted.
        Does nothing on sucesss, raises and error on failure."""
        if self.marate is None and self.propensity_function is None:
            raise ReactionError("You must specify either a mass-action rate or a propensity function")
        if self.marate is None and self.propensity_function is not None:
            if self.propensity_function.strip() == '':
                raise ReactionError("Propensity Function cannot be empty string.")
        if len(self.reactants) == 0 and len(self.products) == 0:
            raise ReactionError("You must have a non-zero number of reactants or products.")
            
    def create_mass_action(self, *args, **kwargs):
        """
        Initializes the mass action propensity function given
        self.reactants and a single parameter value.
        """
        from gillespy2.core import log
        log.warning(
            """
            Reaction.create_mass_action has been deprecated.  Future releases of GillesPy2 may
            not support this feature.
            """
        )

        self.__create_mass_action(*args, **kwargs)





# class Reaction(SortableObject, Jsonify):
#     """
#     Models a single reaction. A reaction has its own dicts of species
#     (reactants and products) and parameters. The reaction's propensity
#     function needs to be evaluable (and result in a non-negative scalar
#     value) in the namespace defined by the union of those dicts.

#     :param name: The name by which the reaction is called (optional).
#     :type name: str

#     :param reactants: The reactants that are consumed in the reaction, with stoichiometry. An
#         example would be {R1 : 1, R2 : 2} if the reaction consumes two of R1 and
#         one of R2, where R1 and R2 are Species objects.
#     :type reactants: dict

#     :param products: The species that are created by the reaction event, with stoichiometry. Same format as reactants.
#     :type products: dict

#     :param propensity_function: The custom propensity function for the reaction. Must be evaluable in the
#         namespace of the reaction using C operations.
#     :type propensity_function: str

#     :param massaction: The switch to use a mass-action reaction. If set to True, a rate value is required.
#     :type massaction: bool

#     :param rate: The rate of the mass-action reaction, take care to note the units.
#     :type rate: float

#     :param annotation: An optional note about the reaction.
#     :type annotation: str

#     Notes
#     ----------
#     For a species that is NOT consumed in the reaction but is part of a mass
#     action reaction, add it as both a reactant and a product.

#     Mass-action reactions must also have a rate term added. Note that the input
#     rate represents the mass-action constant rate independent of volume.

#     if a name is not provided for reactions, the name will be populated by the
#     model based on the order it was added. This could impact seeded simulation 
#     results if the order of addition is not preserved.
#     """

#     def __init__(self, name="", reactants={}, products={}, propensity_function=None, massaction=False, rate=None,
#                  annotation=None):
#         """
#         Initializes the reaction using short-hand notation.
#         """

#         # Metadata
#         if name in (None, ""):
#             self.name = f'rxn{uuid.uuid4()}'.replace('-', '_')
#         else:
#             self.name = name
#         self.annotation = ""

#         # We might use this flag in the future to automatically generate
#         # the propensity function if set to True.
#         if propensity_function is not None:
#             self.massaction = False
#             self.marate = None
#         else:
#             self.massaction = True

#         self.propensity_function = propensity_function
#         self.ode_propensity_function = propensity_function

#         if self.propensity_function is not None and self.massaction:
#             errmsg = ("Reaction {} You cannot set the propensity type to mass-action and simultaneously set a "
#                       "propensity function.").format(self.name)
#             raise ReactionError(errmsg)

#         self.reactants = {}
#         for r in reactants:
#             rtype = type(r).__name__
#             if rtype == 'instance':
#                 self.reactants[r.name] = reactants[r]
#             else:
#                 self.reactants[r] = reactants[r]

#         self.products = {}
#         for p in products:
#             rtype = type(p).__name__
#             if rtype == 'instance':
#                 self.products[p.name] = products[p]
#             else:
#                 self.products[p] = products[p]

#         if self.massaction:
#             self.type = "mass-action"
#             if rate is None:
#                 self.marate = None
#             else:
#                 rtype = type(rate).__name__
#                 if rtype == 'instance':
#                     self.marate = rate.name
#                 else:
#                     self.marate = rate
#                 self.create_mass_action()
#         else:
#             self.type = "customized"

#             def __customPropParser():
#                 pow_func = ast.parse("pow", mode="eval").body

#                 class ExpressionParser(ast.NodeTransformer):
#                     def visit_BinOp(self, node):
#                         node.left = self.visit(node.left)
#                         node.right = self.visit(node.right)
#                         if isinstance(node.op, (ast.BitXor, ast.Pow)):
#                             # ast.Call calls defined function, args include which nodes
#                             # are effected by function call
#                             call = ast.Call(func=pow_func,
#                                             args=[node.left, node.right],
#                                             keywords=[])
#                             # Copy_location copies lineno and coloffset attributes
#                             # from old node to new node. ast.copy_location(new_node,old_node)
#                             call = ast.copy_location(call, node)
#                             # Return changed node
#                             return call
#                         # No modification to node, classes extending NodeTransformer methods
#                         # Always return node or value
#                         else:
#                             return node

#                     def visit_Name(self, node):
#                         # Visits Name nodes, if the name nodes "id" value is 'e', replace with numerical constant
#                         if node.id == 'e':
#                             nameToConstant = ast.copy_location(ast.Num(float(np.e), ctx=node.ctx), node)
#                             return nameToConstant
#                         return node

#                 self.verify()
#                 expr = self.propensity_function
#                 expr = expr.replace('^', '**')
#                 expr = ast.parse(expr, mode='eval')
#                 expr = ExpressionParser().visit(expr)

#                 class ToString(ast.NodeVisitor):
#                     substitutions = {
#                         "ln": "log",
#                     }

#                     def __init__(self):
#                         self.string = ''

#                     def _string_changer(self, addition):
#                         self.string += addition

#                     def visit_BinOp(self, node):
#                         self._string_changer('(')
#                         self.visit(node.left)
#                         self.visit(node.op)
#                         self.visit(node.right)
#                         self._string_changer(')')

#                     def visit_Name(self, node):
#                         self._string_changer(node.id)
#                         self.generic_visit(node)

#                     def visit_Num(self, node):
#                         self._string_changer(str(node.n))
#                         self.generic_visit(node)

#                     def visit_Call(self, node):
#                         func_name = ToString.substitutions.get(node.func.id) \
#                             if node.func.id in ToString.substitutions \
#                             else node.func.id
#                         self._string_changer(func_name + '(')
#                         counter = 0
#                         for arg in node.args:
#                             if counter > 0:
#                                 self._string_changer(',')
#                             self.visit(arg)
#                             counter += 1
#                         self._string_changer(')')

#                     def visit_Add(self, node):
#                         self._string_changer('+')
#                         self.generic_visit(node)

#                     def visit_Div(self, node):
#                         self._string_changer('/')
#                         self.generic_visit(node)

#                     def visit_Mult(self, node):
#                         self._string_changer('*')
#                         self.generic_visit(node)

#                     def visit_UnaryOp(self, node):
#                         self._string_changer('(')
#                         self.visit_Usub(node)
#                         self._string_changer(')')

#                     def visit_Sub(self, node):
#                         self._string_changer('-')
#                         self.generic_visit(node)

#                     def visit_Usub(self, node):
#                         self._string_changer('-')
#                         self.generic_visit(node)

#                 newFunc = ToString()
#                 newFunc.visit(expr)
#                 return newFunc.string

#             self.propensity_function = __customPropParser()
#             self.ode_propensity_function = self.propensity_function

#     def __str__(self):
#         self.verify()
#         print_string = self.name
#         if len(self.reactants):
#             print_string += '\n\tReactants'
#             for r, stoich in self.reactants.items():
#                 try:
#                     if isinstance(r, str):
#                         print_string += '\n\t\t' + r + ': ' + str(stoich)
#                     else:
#                         print_string += '\n\t\t' + r.name + ': ' + str(stoich)
#                 except Exception as e:
#                     print_string += '\n\t\t' + r + ': ' + 'INVALID - ' + str(e)
#         if len(self.products):
#             print_string += '\n\tProducts'
#             for p, stoich in self.products.items():
#                 try:
#                     if isinstance(p, str):
#                         print_string += '\n\t\t' + p + ': ' + str(stoich)
#                     else:
#                         print_string += '\n\t\t' + p.name + ': ' + str(stoich)
#                 except Exception as e:
#                     print_string += '\n\t\t' + p + ': ' + 'INVALID - ' + str(e)
#         print_string += '\n\tPropensity Function: ' + self.propensity_function
#         return print_string

#     def verify(self):
#         """ Check if the reaction is properly formatted.
#         Does nothing on sucesss, raises and error on failure."""
#         if self.marate is None and self.propensity_function is None:
#             raise ReactionError("You must specify either a mass-action rate or a propensity function")
#         if self.marate is None and self.propensity_function is not None:
#             if self.propensity_function.strip() == '':
#                 raise ReactionError("Propensity Function cannot be empty string.")
#         if len(self.reactants) == 0 and len(self.products) == 0:
#             raise ReactionError("You must have a non-zero number of reactants or products.")

#     def create_mass_action(self):
#         """
#         Initializes the mass action propensity function given
#         self.reactants and a single parameter value.
#         """
#         # We support zeroth, first and second order propensities only.
#         # There is no theoretical justification for higher order propensities.
#         # Users can still create such propensities if they really want to,
#         # but should then use a custom propensity.
#         total_stoch = 0
#         for r in sorted(self.reactants):
#             total_stoch += self.reactants[r]
#         if total_stoch > 2:
#             raise ReactionError("Reaction: A mass-action reaction cannot involve more than two of one species or one "
#                                 "of two species. To declare a custom propensity, replace 'rate' with "
#                                 "'propensity_function'.")

#         # Case EmptySet -> Y

#         if isinstance(self.marate, str):
#             propensity_function = self.marate
#             ode_propensity_function = self.marate
#         else:
#             propensity_function = self.marate.name
#             ode_propensity_function = self.marate.name

#         # There are only three ways to get 'total_stoch==2':
#         for r in sorted(self.reactants):
#             if isinstance(r, str):
#                 rname = r
#             else:
#                 rname = r.name
#             # Case 1: 2X -> Y
#             if self.reactants[r] == 2:
#                 propensity_function = "0.5 * {0} * {1} * ({1} - 1) / vol".format(propensity_function, rname)
#                 ode_propensity_function += '*' + rname + '*' + rname
#             else:
#                 # Case 3: X1, X2 -> Y;
#                 propensity_function += "*" + rname
#                 ode_propensity_function += '*' + rname

#         # Set the volume dependency based on order.
#         order = len(self.reactants)
#         if order == 2:
#             propensity_function += "/vol"
#         elif order == 0:
#             propensity_function += "*vol"

#         self.propensity_function = propensity_function
#         self.ode_propensity_function = ode_propensity_function

#     def setType(self, rxntype):
#         """
#         Sets reaction type to either "mass-action" or "customized"

#         :param rxntype: Either "mass-action" or "customized"
#         :type rxntype: str
#         """
#         if rxntype.lower() not in {'mass-action', 'customized'}:
#             raise ReactionError("Invalid reaction type.")
#         self.type = rxntype.lower()

#         self.massaction = False if self.type == 'customized' else True

#     def addReactant(self, S, stoichiometry):
#         """
#         Adds a reactant to the reaction (species that is consumed)

#         :param S: Reactant to add to this reaction.
#         :type S: gillespy2.Species
#         :param stoichiometry: The stoichiometry of the given reactant.
#         :type stoichiometry: int
#         """
#         if stoichiometry <= 0:
#             raise ReactionError("Reaction Stoichiometry must be a \
#                                     positive integer.")
#         self.reactants[S.name] = stoichiometry

#     def addProduct(self, S, stoichiometry):
#         """
#         Adds a product to the reaction (species that is created)

#         :param S: Product to add to this reaction.
#         :type S: gillespy2.Species
#         :param stoichiometry: The stoichiometry of the given product.
#         :type stoichiometry: int
#         """
#         self.products[S.name] = stoichiometry

#     def Annotate(self, annotation):
#         """
#         Adds a note to the reaction

#         :param annotation: An optional note about the reaction.
#         :type annotation: str
#         """
#         self.annotation = annotation

#     def sanitized_propensity_function(self, species_mappings, parameter_mappings, ode=False):
#         names = sorted(list(species_mappings.keys()) + list(parameter_mappings.keys()), key=lambda x: len(x),
#                        reverse=True)
#         replacements = [parameter_mappings[name] if name in parameter_mappings else species_mappings[name]
#                         for name in names]
#         sanitized_propensity = self.ode_propensity_function if ode else self.propensity_function
#         for id, name in enumerate(names):
#             sanitized_propensity = sanitized_propensity.replace(name, "{" + str(id) + "}")
#         return sanitized_propensity.format(*replacements)

#     def to_dict(self):
#         temp = vars(self).copy()

#         # This to_dict overload is needed because, while Species is hashable (thanks to its inheretence),
#         # objects are not valid key values in the JSON spec. To fix this, we set the Species equal to some key 'key', 
#         # and it's value equal to some key 'value'.

#         temp["products"] = list({ "key": k, "value": v} for k, v in self.products.items() )
#         temp["reactants"] = list( {"key": k, "value": v} for k, v in self.reactants.items() )

#         return temp

#     @classmethod
#     def from_json(cls, json_object):
#         new = Reaction.__new__(Reaction)
#         new.__dict__ = json_object

#         # Same as in to_dict(), but we need to reverse it back into its original representation.
#         new.products = { x["key"]: x["value"] for x in json_object["products"] }
#         new.reactants = { x["key"]: x["value"] for x in json_object["reactants"] }

#         return new