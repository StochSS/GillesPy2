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
    def __init__(self, statement: str):
        """
        The provided statement is expected to be a valid Python expression.

        :param statement: Python expression to parse and convert.
        :type statement: str

        :raises SyntaxError: When the provided statement is not a valid Python expression.
        """
        self.statement = statement
        self.language = "python"
        self.expression = Expression.parse_python(statement)
        self.blacklist: "Union[list[ast.operator], None]" = None
        self.namespace: "Union[dict[str, any], None]" = None

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
            print(node.op)
            print(dir(node.op))
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

    def with_blacklist(self, blacklist: "list[ast.operator]" = None) -> "Expression":
        """
        Override operator handling behavior when converting or validating the expression.

        :param blacklist: List of operators which are not allowed.
        :type blacklist: list[ast.operator]

        :returns: List of operators which, based on the given parameters, are not valid.
        An empty list indicates that no invalid operators were found.
        """
        self.blacklist = blacklist
        return self

    def with_namespace(self, namespace: "dict[str, str]") -> "Expression":
        """
        Validates the expression string against the provided namespace.
        Any identifiers present in the expression which are not listed in the namespace
        will cause the expression to be flagged as an invalid namespace.

        :param namespace: A dictionary containing the namespace mappings for the expression.
        The keys of the dict are expected to be the "only" valid identifiers.
        The values of the namespace are what the keys map to during sanitization, if used.
        :type namespace: dict[str, str]

        :returns: List of invalid namespace identifiers found in the expression.
        """
        self.namespace = namespace
        return self

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
