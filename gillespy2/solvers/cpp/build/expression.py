# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import ast
from typing import Union, Optional


class Expression:
    """
    Accepts an expression string to validate and convert.
    Allows for pre-flight syntax and namespace validations,
    as well as converting between Python and C++ expressions.
    """
    def __init__(self,
                 blacklist: "Union[list[str], dict[ast.operator, str]]" = None,
                 namespace: "dict[str, any]" = None,
                 sanitize=False):
        """
        Object for managing context to validate Python expressions.
        Expressions can be passed and validated, which are validated for syntax, namespace, and other conditions.
        The later provided statements are expected to be valid Python expressions.

        :param blacklist: List of operators which are not allowed in the following expressions.
        Note that this will be "forwarded" to all following expressions.
        Ideally, one should define the "universal" blacklist in the constructor,
        using the `Expression#with_blacklist` method for more granular validations.
        :type blacklist: list[str]

        :param namespace: Dictionary mapping allowed bare identifiers to their sanitized equivalents.
        Any bare identifiers not listed as a namespace key will trigger a failed validation.
        :type  namespace: dict[str, any]

        :param sanitize: Whether or not to substitute namespace names during conversion.
        Any valid names found as namespace keys will automatically be converted to the
        corresponding namespace values in the namespace dict when getexpr_* methods are called.
        :type sanitize: bool
        """
        if blacklist is None:
            blacklist = dict({})
        elif not isinstance(blacklist, dict):
            blacklist = {ast_op: op for op, ast_op in Expression.map_operator(blacklist)}
        if namespace is None:
            namespace = {}
        self.blacklist = blacklist
        self.namespace = namespace
        self.sanitize = sanitize

    class ValidationVisitor(ast.NodeTransformer):
        def __init__(self,
                     namespace: "dict[str, str]" = None,
                     blacklist: "dict[ast.operator, str]" = None,
                     sanitize=False):
            self.namespace = dict({}) if namespace is None else namespace
            self.blacklist = dict({}) if blacklist is None else blacklist
            self.sanitize = sanitize
            self.invalid_names = set()
            self.invalid_operators = set()

        def check_blacklist(self, operator):
            operator = type(operator)
            if operator in self.blacklist:
                self.invalid_operators.add(str(self.blacklist.get(operator)))

        def visit_Name(self, node: "ast.Name"):
            if node.id not in self.namespace:
                self.invalid_names.add(node.id)
            elif self.sanitize:
                node.id = self.namespace.get(node.id)
            self.generic_visit(node)
            return node

        def visit_Call(self, node: "ast.Call"):
            if node.func.id not in self.namespace:
                self.invalid_names.add(node.func.id)
            elif self.sanitize:
                node.func.id = self.namespace.get(node.func.id)
            self.generic_visit(node)
            return node

        def visit_BinOp(self, node: "ast.BinOp"):
            self.check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_UnaryOp(self, node: "ast.UnaryOp"):
            self.check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_BoolOp(self, node: "ast.BoolOp"):
            self.check_blacklist(node.op)
            self.generic_visit(node)
            return node

        def visit_Compare(self, node: "ast.Compare"):
            for op in node.ops:
                self.check_blacklist(op)
            self.generic_visit(node)
            return node

        def visit_Assign(self, node: "ast.Assign"):
            self.check_blacklist(ast.Assign())
            self.generic_visit(node)
            return node

    operator_map = {
        # Basic math operators
        "+": ast.Add,
        "-": ast.Sub,
        "*": ast.Mult,
        "/": ast.Div,
        "**": ast.Pow,
        "//": ast.FloorDiv,
        "%": ast.Mod,
        "@": ast.MatMult,

        # Variable operators
        "=": ast.Assign,
        ":=": ast.Assign,

        # Boolean operators
        "and": ast.And,
        "or": ast.Or,

        # Bitwise operators (^ gets substituted for ** in Python and pow() in C++)
        "^": ast.BitXor,
        "<<": ast.LShift,
        ">>": ast.RShift,
        "|": ast.BitOr,
        "&": ast.BitAnd,

        # Comparison operators
        "==": ast.Eq,
        "!=": ast.NotEq,
        "!": ast.Not,
        ">": ast.Gt,
        ">=": ast.GtE,
        "<": ast.Lt,
        "<=": ast.LtE,
    }

    @classmethod
    def map_operator(cls, operator: "Union[str, list[str]]"):
        if isinstance(operator, list):
            for op in operator:
                yield from Expression.map_operator(op)
        else:
            # Base case: operator is a single string.
            if operator in Expression.operator_map:
                yield operator, Expression.operator_map.get(operator)
            elif operator in Expression.operator_map.values():
                # Yield the operator directly if there is no need to map it.
                yield operator

    def with_blacklist(self, blacklist: "list[str]" = None) -> "Expression":
        """
        Create a new duplicate of the current expression, with a different operator blacklist.
        Overrides operator handling behavior when converting or validating the expression.

        :param blacklist: List of operators which are not allowed.
        :type blacklist: list[ast.operator]

        :returns: List of operators which, based on the given parameters, are not valid.
        An empty list indicates that no invalid operators were found.
        """
        if blacklist is None:
            blacklist = self.blacklist.copy()
        return Expression(blacklist=blacklist, namespace=self.namespace)

    def with_namespace(self, namespace: "dict[str, any]" = None) -> "Expression":
        """
        Create a new duplicate of the current expression, with a different namespace.
        Any identifiers present in the expression which are not listed in the namespace
        will cause the expression to be flagged as an invalid namespace during validation.

        :param namespace: A dictionary containing the namespace mappings for the expression.
        The keys of the dict are expected to be the "only" valid identifiers.
        The values of the namespace are what the keys map to during sanitization, if used.
        :type namespace: dict[str, str]

        :returns: New expression containing the given namespace.
        The returned expression is a *copy* of the current expression.
        :rtype: Expression
        """
        if namespace is None:
            namespace = self.namespace.copy()
        return Expression(blacklist=self.blacklist.copy(), namespace=namespace)

    def validate(self, statement: "str") -> "ExpressionResults":
        """
        Using the information provided so far, ensure that the given Python expression is valid.
        The Python expression is parsed, raising a SyntaxError if it is an invalid Python expression.
        The expression is then checked against the given properties, such as namespace and operator blacklist.
        Additionally, the expression is rejected if it is not a single rvalue expression.

        :raises SyntaxError: The statement is not a valid Python expression.
        :returns: True if the statement is valid, otherwise returns false.
        """
        expr = ast.parse(statement)
        validator = Expression.ValidationVisitor(self.namespace, self.blacklist, self.sanitize)
        validator.visit(expr)

        return ExpressionResults(invalid_names=validator.invalid_names, invalid_operators=validator.invalid_operators)

    def __get_expr(self, statement: "str", converter: "ExpressionConverter") -> "Optional[str]":
        validator = Expression.ValidationVisitor(self.namespace, self.blacklist, self.sanitize)
        validator.visit(converter.tree)

        failures_found = []

        if validator.invalid_operators:
            base_msg = "Blacklisted operator"
            base_msg = f"{base_msg}s" if len(validator.invalid_operators) > 1 else base_msg
            failures_found.append(f"{base_msg}: {','.join(validator.invalid_operators)}")

        if validator.invalid_names:
            base_msg = "Cannot resolve species name"
            base_msg = f"{base_msg}s" if len(validator.invalid_names) > 1 else base_msg
            failures_found.append(f"{base_msg}: {','.join(validator.invalid_names)}")

        if len(failures_found) > 0:
            raise SyntaxError(f"Invalid GillesPy2 expression \"{statement}\"\n"
                              + "\n".join([f"* {msg}" for msg in failures_found]))

        return converter.get_str()

    def getexpr_python(self, statement: "str") -> "Optional[str]":
        """
        Converts the expression object into a Python expression string.
        Raises a SyntaxError if conversion to a Python string is impossible.

        :raises: SyntaxError
        :returns: Python expression string, if valid. Returns None if validation fails.
        """
        expr = ast.parse(statement)
        return self.__get_expr(statement, PythonConverter(expr))

    def getexpr_cpp(self, statement: "str") -> "Optional[str]":
        """
        Converts the expression object into a C++ expression string.
        Raises a SyntaxError if conversion to a C++ string is impossible.

        :raises: SyntaxError
        :returns: C++ expression string
        """
        statement = ExpressionConverter.convert_str(statement)
        expr = ast.parse(statement)
        return self.__get_expr(statement, CppConverter(expr))


class ExpressionResults:
    """
    Container struct for returning the results of expression validation.
    Any expression items which indicate an invalid expression are listed on an ExpressionResults instance.
    Empty lists indicate that the expression is valid.
    """
    def __init__(self, invalid_names: "set[str]" = None, invalid_operators: "set[str]" = None, is_valid=True):
        """
        Container struct for returning the results of expression validation.

        :param invalid_names: List of expression identifiers which were not valid in the given namespace.
        :type invalid_names: list[str]

        :param invalid_operators: List of blacklisted operators which were present in the expression.
        :type invalid_operators: list[str]

        :param is_valid: Override value for the `is_valid` property.
        If not set, then the validity of the expression is inferred by the `invalid_*` lists provided.
        """
        self.invalid_names = invalid_names
        self.invalid_operators = invalid_operators
        self.is_valid = is_valid and (not invalid_names and not invalid_operators)


class ExpressionConverter(ast.NodeVisitor):
    def __init__(self, tree: "ast.AST"):
        self.tree = tree
        self.expression = []

    @classmethod
    def convert_str(cls, expression: "str"):
        return expression.replace("^", "**")

    def parse_operator(self, operator: "str"):
        expr = f"({self.expression.pop()}{operator}{self.expression.pop()})"
        self.expression.append(expr)

    def parse_logical(self, operator: "str"):
        expr = f"{self.expression.pop()} {operator} {self.expression.pop()}"
        self.expression.append(expr)

    def parse_comparison(self, comparator: "str"):
        expr = f"{self.expression.pop()} {comparator} {self.expression.pop()}"
        self.expression.append(expr)

    def visit_Name(self, node: "ast.Name"):
        self.expression.append(node.id)
        self.generic_visit(node)

    def visit_Constant(self, node: "ast.Constant"):
        self.expression.append(str(node.value))
        self.generic_visit(node)

    ###########################################################################
    ### The below methods are deprecated as of Python 3.8.                  ###
    ### They are included for compatibility with 3.7 and earlier.           ###
    ### If <=3.7 becomes unsupported, please remove these visitor methods.  ###
    def visit_Num(self, node: "ast.Num"):
        self.expression.append(str(node.n))
        self.generic_visit(node)

    def visit_Str(self, node: "ast.Str"):
        self.expression.append(str(node.s))
        self.generic_visit(node)

    def visit_Bytes(self, node: "ast.Bytes"):
        self.expression.append(str(node.s))
        self.generic_visit(node)

    def visit_NameConstant(self, node: "ast.NameConstant"):
        self.expression.append(str(node.value))
        self.generic_visit(node)

    def visit_Ellipsis(self, node: "ast.Ellipsis"):
        self.expression.append(str(node))
        self.generic_visit(node)
    ### End of deprecated functions (deprecated as of Python 3.8)           ###
    ###########################################################################

    def visit_BinOp(self, node: "ast.BinOp"):
        # Right node is visited first.
        # By visiting the left node last, the most recently appended token is always the left-hand token.
        # This allows us to always append when adding to the expression, and always pop when processing it.
        self.visit(node.right)
        self.visit(node.left)
        self.visit(node.op)

    def visit_Call(self, node: "ast.Call"):
        arg_list = []
        for arg in node.args:
            self.visit(arg)
            arg_list.append(self.expression.pop())
        arg_list = ",".join(arg_list)
        expr = f"{node.func.id}({arg_list})"
        self.expression.append(expr)

    def visit_BoolOp(self, node: "ast.BoolOp"):
        # Base converter class assumes that "And" and "Or" operations are defined by inheriting class.
        # Implement visit_And() and visit_Or() to define behavior.
        for operand in reversed(node.values):
            self.visit(operand)
        # Process n-1 operations; for n operands, there are n-1 operations.
        # Example: x && y || z -> 3 operands, 2 operations
        for op_i in range(len(node.values) - 1):
            self.visit(node.op)

    def visit_Add(self, node: "ast.Add"):
        self.generic_visit(node)
        self.parse_operator("+")

    def visit_Sub(self, node: "ast.Sub"):
        self.generic_visit(node)
        self.parse_operator("-")

    def visit_Mult(self, node: "ast.Mult"):
        self.generic_visit(node)
        self.parse_operator("*")

    def visit_Div(self, node: "ast.Div"):
        self.generic_visit(node)
        self.parse_operator("/")

    def visit_Pow(self, node: "ast.Pow"):
        self.generic_visit(node)
        self.parse_operator("**")

    def visit_Compare(self, node: "ast.Compare"):
        for comparator in node.comparators:
            self.visit(comparator)
        self.visit(node.left)
        for operator in node.ops:
            self.visit(operator)

    def visit_Eq(self, node: "ast.Eq"):
        self.generic_visit(node)
        self.parse_comparison("==")

    def visit_NotEq(self, node: "ast.NotEq"):
        self.generic_visit(node)
        self.parse_comparison("!=")

    def visit_Lt(self, node: "ast.Lt"):
        self.generic_visit(node)
        self.parse_comparison("<")

    def visit_LtE(self, node: "ast.LtE"):
        self.generic_visit(node)
        self.parse_comparison("<=")

    def visit_Gt(self, node: "ast.Gt"):
        self.generic_visit(node)
        self.parse_comparison(">")

    def visit_GtE(self, node: "ast.GtE"):
        self.generic_visit(node)
        self.parse_comparison(">=")

    def visit_UnaryOp(self, node: "ast.UnaryOp"):
        self.visit(node.operand)
        self.visit(node.op)

    def visit_USub(self, node: "ast.USub"):
        self.expression.append(f"-{self.expression.pop()}")

    def _get_str(self, expr: "ast.AST"):
        self.visit(expr)
        return "".join(self.expression)

    def get_str(self) -> "str":
        return self._get_str(self.tree)


class PythonConverter(ExpressionConverter):
    def visit_And(self, node: "ast.And"):
        self.parse_logical("and")

    def visit_Or(self, node: "ast.Or"):
        self.parse_logical("or")


class CppConverter(ExpressionConverter):
    class CppExpressionTransformer(ast.NodeTransformer):
        def visit_BinOp(self, node: "ast.BinOp"):
            self.generic_visit(node)
            if isinstance(node.op, ast.Pow):
                node = ast.copy_location(ast.Call(
                    func=ast.Name(id='pow', ctx=ast.Load()),
                    args=[node.left, node.right],
                    keywords=[]
                ), node)
            return node

    def visit_And(self, node: "ast.And"):
        self.generic_visit(node)
        self.parse_logical("&&")

    def visit_Or(self, node: "ast.Or"):
        self.generic_visit(node)
        self.parse_logical("||")

    def get_str(self) -> "str":
        expr = CppConverter.CppExpressionTransformer().visit(self.tree)
        return super()._get_str(expr)
