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

import unittest
import numpy
import io
import os
from example_models import create_dimerization
import subprocess
import tempfile
from gillespy2.solvers.cpp.c_decoder import BasicSimDecoder
from gillespy2.solvers.cpp import SSACSolver, ODECSolver, TauLeapingCSolver
from gillespy2.solvers.cpp import TauHybridCSolver
from gillespy2.solvers.cpp.build.expression import Expression, ExpressionConverter
from gillespy2.solvers.cpp.build.template_gen import SanitizedModel


class ExpressionTestCase:
    """
    Each test expression consists of a dict of argument names, an expression, and a list of
    values to be passed as arguments to the given expression.
    """
    def __init__(self, args: "dict[str, str]", expression: "str", values: "list[list[float]]"):
        self.args = args
        self.expression = expression
        self.values = values


class TestCSolvers(unittest.TestCase):
    """
    Run tests to ensure the setup process is valid for each C++ solver.
    Does not actually run any simulations, just validates their construction.
    Effectively an "integration test" for each C++ solver, making sure they
      play nice with the build engine and template gen.
    """

    test_model = create_dimerization()
    base_solvers = [
        SSACSolver,
        ODECSolver,
        TauLeapingCSolver,
        TauHybridCSolver,
    ]
    solvers = {
        SSACSolver.target: SSACSolver(model=test_model),
        ODECSolver.target: ODECSolver(model=test_model),
        TauLeapingCSolver.target: TauLeapingCSolver(model=test_model),
        TauHybridCSolver.target: TauHybridCSolver(model=test_model)
    }
    solvers_variable = {
        SSACSolver.target: SSACSolver(model=test_model, variable=True),
        ODECSolver.target: ODECSolver(model=test_model, variable=True),
        TauLeapingCSolver.target: TauLeapingCSolver(model=test_model, variable=True),
        TauHybridCSolver.target: TauHybridCSolver(model=test_model, variable=True),
    }
    expressions = [
        # Each test expression consists of a dict of args, an expression string, and a list of arg values.
        # Asserts that single operations work.
        ExpressionTestCase({"x": "x"}, "x*2", [
            [0.0], [1.0], [-1.0], [9.999], [-9.999],
        ]),
        # Asserts that order of operations is being evaluated properly.
        ExpressionTestCase({"x": "x"}, "x*2 + x/2 - (x*3)^2 + x/3^2", [
            [0.0], [1.0], [-1.0], [3.333], [-3.333], [9.8765], [-9.8765]
        ]),
        # Asserts that order of operations is evaluated properly with multiple variables.
        ExpressionTestCase({"x": "x", "y": "y"}, "(x-1)*y^2+x", [
            [1.0, 2.4], [5.1, 0.0], [5.1, 1.0], [5.1, -1.0], [9.8765, -1.0], [-1.0, 9.8765],
        ]),
        # Asserts complex order of operations with a large number of variables.
        ExpressionTestCase({"x": "x", "y": "y", "z": "z"}, "(x^2/y^2/z^2)/x^2/y^2/z^2**1/x**1/y**1/z", [
            [5.1, 0.1, 2.0], [0.1, 5.1, 2.0], [2.0, 0.1, 5.1], [2.0, 5.1, 0.1],
        ]),
        # Known, builtin math expression functions work.
        ExpressionTestCase({"x": "x"}, "abs(x)", [
            [100.0], [100], [-100.0], [-100], [0],
        ]),
    ]
    comparisons = [
        # Asserts that single comparison expressions work.
        ExpressionTestCase({"x": "x"}, "x > 0", [
            [100], [0], [0.001], [-1],
        ]),
        ExpressionTestCase({"x": "x", "y": "y"}, "x > y", [
            [100, 99], [99, 100], [-10, 10], [10, -10],
            [0.001, 0.0], [0.0, 0.001], [-99.999, -99.998]
        ]),
        # Asserts that single boolean operators work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < x", [
            [100, 99], [99, 100], [0, -100], [-0.001, -99.0], [0, 0.001], [-0.001, 0]
        ]),
        # Asserts that nested boolean operators work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < 10 and x > y", [
            [100, 9], [0.01, 0.00], [100, 200], [0.01, 0.02],
            [0, 0], [-0.01, -0.02], [-0.01, 0],
        ]),
        # Asserts that both && and || work.
        ExpressionTestCase({"x": "x", "y": "y"}, "x > 0 and y < 10 or y > 100", [
            [10, 9], [0.01, 9.99], [0, 10], [-1.0, -1.0],
        ]),
        # Asserts that nested boolean operators properly respect order of operations.
        ExpressionTestCase({"x": "x", "y": "y", "z": "z"}, "x^2>x and y<y^2 or z^2!=z^3 and y!=z", [
            [1.0, 1.0, 1.0], [99.9, 99.9, 100.0],
            [0.0, -1.0, 99.9], [-1.0, -1.0, 0.00],
        ]),
    ]

    def test_solver_build(self):
        """
        Build each solver and ensure that they build properly.
        """
        # Test builds for non-variable solvers
        for solver_name, solver in self.solvers.items():
            with self.subTest(solver=solver_name):
                exe = solver._build(model=self.test_model,
                                    simulation_name=solver_name,
                                    variable=False,
                                    debug=False)

                self.assertTrue(os.path.isfile(exe),
                                "Built simulation output could not be found or is a directory.")
                self.assertTrue(os.access(exe, os.X_OK),
                                "Built simulation binaries are not executable.")
                solver.build_engine.clean()
                self.assertFalse(os.path.exists(exe),
                                 "Simulation output not cleaned up after call to .clean().")

        # Test builds for variable solvers
        for solver_name, solver in self.solvers_variable.items():
            with self.subTest(solver=solver_name):
                exe = solver._build(model=self.test_model,
                                    simulation_name=solver_name,
                                    variable=True,
                                    debug=False)
                self.assertTrue(os.path.isfile(exe),
                                "Built simulation output could not be found or is a directory.")
                self.assertTrue(os.access(exe, os.X_OK),
                                "Built simulation binaries are not executable.")
                solver.build_engine.clean()
                self.assertFalse(os.path.exists(exe),
                                 "Simulation output not cleaned up after call to .clean().")

    def test_solver_precompile(self):
        """
        Ensure that a pre-built C++ solver actually builds the simulation before running.
        """
        for solver in self.base_solvers:
            with self.subTest(solver=solver.name):
                base_solver = solver(model=self.test_model)
                self.assertIsNotNone(base_solver.build_engine, 
                                     "Build engine not created at solver's construction")
                self.assertIsNotNone(base_solver.build_engine.get_executable_path(),
                                     "Build engine has no associated executable")
                self.assertTrue(os.access(base_solver.build_engine.get_executable_path(), os.X_OK),
                                "Solver executable invalid or missing at solver's construction")

    def test_solver_expressions(self):
        """
        Ensure that expression conversions to C++ result in (roughly) equivalent values as Python.
        """
        tmpdir = tempfile.mkdtemp()
        src_path = os.path.join(os.path.dirname(__file__), "assets", "evaluate.c")
        exe_path = os.path.join(tmpdir, "test")

        def build(expr_args: "list[str]", expr_str: "str", use_bool=False):
            args = ["gcc", "-o", exe_path, src_path, "-lm"]
            expr_num = str(len(expr_args))
            expr_args = ",".join(expr_args)
            args.append(f"-DEXP{expr_num}({expr_args})=({expr_str})")
            if use_bool:
                args.append("-DUSE_BOOLEAN")
            subprocess.check_call(args)

        def run(args: "list[str]") -> str:
            args.insert(0, exe_path)
            stdout = subprocess.check_output(args)
            return stdout.decode("ascii")

        def test_expressions(expressions: "list[ExpressionTestCase]", use_bool=False):
            for entry in expressions:
                expression = ExpressionConverter.convert_str(entry.expression)
                expr = Expression(namespace={
                    **SanitizedModel.reserved_names,
                    **SanitizedModel.function_map,
                    **entry.args,
                })
                cpp_expr = expr.getexpr_cpp(expression)
                with self.subTest(msg="Evaluating converted C expressions",
                                  expression=entry.expression,
                                  c_expression=cpp_expr):
                    py_args = ",".join(entry.args.keys())
                    py_func = eval(f"lambda {py_args}: {expression}")

                    for value_set in entry.values:
                        value_str = [str(val) for val in value_set]
                        with self.subTest(values=",".join(value_str)):
                            expect = py_func(*value_set)
                            build(list(entry.args.values()), cpp_expr, use_bool)
                            if use_bool:
                                result_cpp = bool(int(run(value_str)))
                                self.assertTrue(expect == result_cpp)
                            else:
                                result_cpp = float(run(value_str))
                                self.assertAlmostEqual(expect, result_cpp, places=3)

        try:
            # Test expressions which return a float value
            test_expressions(self.expressions, False)
            # Test boolean and comparator expressions
            test_expressions(self.comparisons, True)

        finally:
            import shutil
            shutil.rmtree(tmpdir)
