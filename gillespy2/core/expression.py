import ast
from gillespy2 import log
from typing import Union


class ExpressionLanguage:
    """Enum class representing the supported languages an Expression object can convert to."""
    PYTHON = "python"
    CPP = "cpp"


class Expression:
    """
    Accepts an expression string to validate and convert.
    Allows for pre-flight syntax and namespace validations,
    as well as converting between Python and C++ expressions.
    """
    def __init__(self, blacklist: "list[str]" = None, namespace: "dict[str, any]" = None):
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
        """
        if blacklist is None:
            blacklist = []
        if namespace is None:
            namespace = {}
        self.language = "python"
        self.blacklist = [op for op in Expression.map_operator(blacklist)]
        self.namespace = namespace

    class ValidationVisitor(ast.NodeVisitor):
        def __init__(self, namespace: "dict[str, any]" = None, blacklist: "list[ast.operator]" = None):
            self.namespace = {} if namespace is None else namespace
            self.blacklist = [] if blacklist is None else blacklist
            self.invalid_names = []
            self.invalid_operators = []

        def visit_Name(self, node: "ast.Name"):
            if node.id not in self.namespace:
                self.invalid_names.append(node.id)

        def visit_BinOp(self, node: "ast.BinOp"):
            if type(node.op) in self.blacklist:
                self.invalid_operators.append(str(node.op))
            self.generic_visit(node)

    @classmethod
    def parse_python(cls, statement: str) -> "Union[ast.AST, None]":
        """
        Attempt to parse the given string as a Python expression.

        :returns: If valid, returns a Python AST representing the parsed expression.
        If the statement is not valid, returns None.
        """
        parsed_expression = None
        try:
            parsed_expression = ast.parse(source=statement)
        except SyntaxError as err:
            log.warning(f"Syntax error: {err.msg}")
        finally:
            return parsed_expression

    operator_map = {
        "+": ast.Add,
        "-": ast.Sub,
        "*": ast.Mult,
        "/": ast.Div,
    }

    @classmethod
    def map_operator(cls, operator: "Union[str, list[str]]"):
        if isinstance(operator, list):
            for op in operator:
                yield from Expression.map_operator(op)
        else:
            # Base case: operator is a single string.
            if operator in Expression.operator_map:
                yield Expression.operator_map[operator]
            elif operator in Expression.operator_map.values():
                # Yield the operator directly if there is no need to map it.
                yield operator

    def with_blacklist(self, blacklist: "list[ast.operator]" = None) -> "Expression":
        """
        Create a new duplicate of the current expression, with a different operator blacklist.
        Overrides operator handling behavior when converting or validating the expression.

        :param blacklist: List of operators which are not allowed.
        :type blacklist: list[ast.operator]

        :returns: List of operators which, based on the given parameters, are not valid.
        An empty list indicates that no invalid operators were found.
        """
        if blacklist is None:
            blacklist = self.blacklist
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
            namespace = self.namespace
        return Expression(blacklist=self.blacklist, namespace=namespace)

    def validate(self, statement: "str") -> "bool":
        """
        Using the information provided so far, ensure that the given Python expression is valid.
        The Python expression is parsed, raising a SyntaxError if it is an invalid Python expression.
        The expression is then checked against the given properties, such as namespace and operator blacklist.
        Additionally, the expression is rejected if it is not a single rvalue expression.

        :raises SyntaxError: The statement is not a valid Python expression.
        """
        expr = ast.parse(statement)

        validator = Expression.ValidationVisitor(namespace=self.namespace, blacklist=self.blacklist)
        validator.visit(expr)

        if validator.invalid_operators:
            log.error("Invalid operators")
            return False

        if validator.invalid_names:
            log.error("Invalid names")
            return False

        return True

    def getexpr_python(self, sanitize=False) -> str:
        """
        Converts the expression object into a Python expression string.
        Raises a SyntaxError if conversion to a Python string is impossible.

        :raises: SyntaxError
        :returns: Python expression string
        """
        raise NotImplementedError("Converting expression to Python string has not yet been implemented.")

    def getexpr_cpp(self, sanitize=False) -> str:
        """
        Converts the expression object into a C++ expression string.
        Raises a SyntaxError if conversion to a C++ string is impossible.

        :raises: SyntaxError
        :returns: C++ expression string
        """
        raise NotImplementedError("Converting expression to C++ string has not yet been implemented.")
