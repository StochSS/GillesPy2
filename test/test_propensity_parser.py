import unittest
from gillespy2.solvers.utilities.solverutils import dependency_grapher
from gillespy2.core import Reaction

r1 = Reaction(name='r1', propensity_function="5*x^2+e*b+6")
r2 = Reaction(name='r2', propensity_function="5*x**2+e*b+6")

r3 = Reaction(name='r3', propensity_function="1*alpha/2+5^beta")
r4 = Reaction(name='r4', propensity_function="1*alpha/2+5**beta")

r5 = Reaction(name='r5', propensity_function="2.78*x+3^(4*x)")
r6 = Reaction(name='r6', propensity_function="2.78*x+3**(4*x)")

r7 = Reaction(name='r7', propensity_function="(alpha/beta + delta**gamma)/(atlas-zeta)")
r8 = Reaction(name='r8', propensity_function="(alpha/beta + delta^gamma)/(atlas-zeta)")

r9 = Reaction(name='r9', propensity_function="-5*-x^2")
r10 = Reaction(name='r10', propensity_function="-5*-x**2")


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
