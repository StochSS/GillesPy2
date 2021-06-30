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

import unittest
from gillespy2.solvers.utilities.solverutils import dependency_grapher
from gillespy2.core import Reaction, Species

s = Species(name='s', initial_value=0)
r1 = Reaction(name='r1', reactants={s:1}, propensity_function="5*x^2+e*b+6")
r2 = Reaction(name='r2', reactants={s:1}, propensity_function="5*x**2+e*b+6")

r3 = Reaction(name='r3', reactants={s:1}, propensity_function="1*alpha/2+5^beta")
r4 = Reaction(name='r4', reactants={s:1}, propensity_function="1*alpha/2+5**beta")

r5 = Reaction(name='r5', reactants={s:1}, propensity_function="2.78*x+3^(4*x)")
r6 = Reaction(name='r6', reactants={s:1}, propensity_function="2.78*x+3**(4*x)")

r7 = Reaction(name='r7', reactants={s:1}, propensity_function="(alpha/beta + delta**gamma)/(atlas-zeta)")
r8 = Reaction(name='r8', reactants={s:1}, propensity_function="(alpha/beta + delta^gamma)/(atlas-zeta)")

r9 = Reaction(name='r9', reactants={s:1}, propensity_function="-5*-x^2")
r10 = Reaction(name='r10', reactants={s:1}, propensity_function="-5*-x**2")


class TestPropensityFunctions(unittest.TestCase):
    def test_propensity_functions(self):
        self.assertEqual(r1.propensity_function, "(((5*pow(x,2))+(2.718281828459045*b))+6)",
                         msg="Has incorrect expression")
        self.assertEqual(r2.propensity_function, "(((5*pow(x,2))+(2.718281828459045*b))+6)",
                         msg="Has incorrect expression")

        self.assertEqual(r3.propensity_function, "(((1*alpha)/2)+pow(5,beta))", msg="Has incorrect expression")
        self.assertEqual(r4.propensity_function, "(((1*alpha)/2)+pow(5,beta))", msg="Has incorrect expression")

        self.assertEqual(r5.propensity_function, "((2.78*x)+pow(3,(4*x)))", msg="Has incorrect expression")
        self.assertEqual(r6.propensity_function, "((2.78*x)+pow(3,(4*x)))", msg="Has incorrect expression")

        self.assertEqual(r7.propensity_function, "(((alpha/beta)+pow(delta,gamma))/(atlas-zeta))",
                         msg="Has incorrect expression")
        self.assertEqual(r8.propensity_function, "(((alpha/beta)+pow(delta,gamma))/(atlas-zeta))",
                         msg="Has incorrect expression")

        self.assertEqual(r9.propensity_function, "((-5)*(-pow(x,2)))",
                         msg="Has incorrect expression")
        self.assertEqual(r10.propensity_function, "((-5)*(-pow(x,2)))",
                         msg="Has incorrect expression")

    def test_dependency_graphing(self):
        from example_models import ToggleSwitch, MichaelisMenten
        model = ToggleSwitch()
        dependencies = dependency_grapher(model, list(model.listOfReactions.keys()))
        correct_graph = {'cu': {'dependencies': ['cv', 'dv']}, 'cv': {'dependencies': ['cu', 'du']},
                         'du': {'dependencies': ['cu']}, 'dv': {'dependencies': ['cv']}}
        self.assertEqual(correct_graph, dependencies)

        model = MichaelisMenten()
        dependencies = dependency_grapher(model, list(model.listOfReactions.keys()))
        correct_graph = {'r1': {'dependencies': ['r2', 'r3']}, 'r2': {'dependencies': ['r1', 'r3']},
                         'r3': {'dependencies': ['r1', 'r2']}}
        self.assertEqual(correct_graph, dependencies)
