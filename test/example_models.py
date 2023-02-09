# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2023 GillesPy2 developers.

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
import gillespy2
from gillespy2.core import (
    Event,
    Model,
    Species,
    Reaction,
    RateRule,
    Parameter,
    EventTrigger,
    AssignmentRule,
    EventAssignment,
	FunctionDefinition
)

import numpy as np


def create_vilar_oscillator(parameter_values=None):
    model = Model(name="VilarOscillator")
    model.volume = 1

    # Parameters
    alphaA = Parameter(name="alphaA", expression=50)
    alphaA_prime = Parameter(name="alphaA_prime", expression=500)
    alphaR = Parameter(name="alphaR", expression=0.01)
    alphaR_prime = Parameter(name="alphaR_prime", expression=50)
    betaA = Parameter(name="betaA", expression=50)
    betaR = Parameter(name="betaR", expression=5)
    deltaMA = Parameter(name="deltaMA", expression=10)
    deltaMR = Parameter(name="deltaMR", expression=0.5)
    deltaA = Parameter(name="deltaA", expression=1)
    deltaR = Parameter(name="deltaR", expression=0.2)
    gammaA = Parameter(name="gammaA", expression=1)
    gammaR = Parameter(name="gammaR", expression=1)
    gammaC = Parameter(name="gammaC", expression=2)
    thetaA = Parameter(name="thetaA", expression=50)
    thetaR = Parameter(name="thetaR", expression=100)
    
    model.add_parameter([alphaA, alphaA_prime, alphaR, alphaR_prime, betaA, betaR,
                        deltaMA, deltaMR, deltaA, deltaR, gammaA, gammaR, gammaC,
                        thetaA, thetaR])

    # Species
    Da = Species(name="Da", initial_value=1, mode="discrete")
    Da_prime = Species(name="Da_prime", initial_value=0, mode="discrete")
    Ma = Species(name="Ma", initial_value=0, mode="discrete")
    Dr = Species(name="Dr", initial_value=1, mode="discrete")
    Dr_prime = Species(name="Dr_prime", initial_value=0, mode="discrete")
    Mr = Species(name="Mr", initial_value=0, mode="discrete")
    C = Species(name="C", initial_value=0, mode="discrete")
    A = Species(name="A", initial_value=0, mode="discrete")
    R = Species(name="R", initial_value=0, mode="discrete")
    
    model.add_species([Da, Da_prime, Ma, Dr, Dr_prime, Mr, C, A, R])

    # Reactions
    r1 = Reaction(name="r1", reactants={'A': 1, 'R': 1}, products={'C': 1}, rate="gammaC")
    r2 = Reaction(name="r2", reactants={'A': 1}, products={}, rate="deltaA")
    r3 = Reaction(name="r3", reactants={'C': 1}, products={'R': 1}, rate="deltaA")
    r4 = Reaction(name="r4", reactants={'R': 1}, products={}, rate="deltaR")
    r5 = Reaction(name="r5", reactants={'A': 1, 'Da': 1}, products={'Da_prime': 1}, rate="gammaA")
    r6 = Reaction(name="r6", reactants={'Da_prime': 1}, products={'A': 1, 'Da': 1}, rate="thetaA")
    r7 = Reaction(name="r7", reactants={'Da': 1}, products={'Da': 1, 'Ma': 1}, rate="alphaA")
    r8 = Reaction(name="r8", reactants={'Da_prime': 1}, products={'Da_prime': 1, 'Ma': 1}, rate="alphaA_prime")
    r9 = Reaction(name="r9", reactants={'Ma': 1}, products={}, rate="deltaMA")
    r10 = Reaction(name="r10", reactants={'Ma': 1}, products={'A': 1, 'Ma': 1}, rate="betaA")
    r11 = Reaction(name="r11", reactants={'A': 1, 'Dr': 1}, products={'Dr_prime': 1}, rate="gammaR")
    r12 = Reaction(name="r12", reactants={'Dr_prime': 1}, products={'A': 1, 'Dr': 1}, rate="thetaR")
    r13 = Reaction(name="r13", reactants={'Dr': 1}, products={'Dr': 1, 'Mr': 1}, rate="alphaR")
    r14 = Reaction(name="r14", reactants={'Dr_prime': 1}, products={'Dr_prime': 1, 'Mr': 1}, rate="alphaR_prime")
    r15 = Reaction(name="r15", reactants={'Mr': 1}, products={}, rate="deltaMR")
    r16 = Reaction(name="r16", reactants={'Mr': 1}, products={'Mr': 1, 'R': 1}, rate="betaR")
    
    model.add_reaction([r1, r2, r3, r4, r5, r6, r7, r8, r9,
                       r10, r11, r12, r13, r14, r15, r16])

    # Timespan
    model.timespan(np.linspace(0,200,201))
    return model


def create_dimerization(parameter_values=None):
    # First call the Model initializer.
    model = Model(name="Dimerization")

    # Define parameters for the rates of creation and dissociation.
    k_c = Parameter(name='k_c', expression=0.005)
    k_d = Parameter(name='k_d', expression=0.08)
    model.add_parameter([k_c, k_d])

    # Define variables for the molecular species representing M and D.
    m = Species(name='monomer', initial_value=30)
    d = Species(name='dimer',   initial_value=0)
    model.add_species([m, d])

    # Define the reactions representing the process.  In GillesPy2,
    # the list of reactants and products for a Reaction object are each a
    # Python dictionary in which the dictionary keys are Species objects
    # and the values are stoichiometries of the species in the reaction.
    r_creation = Reaction(name="r_creation", rate=k_c,
                                        reactants={m:2}, products={d:1})
    r_dissociation = Reaction(name="r_dissociation", rate=k_d,
                                        reactants={d:1}, products={m:2})
    model.add_reaction([r_creation, r_dissociation])

    # Set the timespan for the simulation.
    model.timespan(np.linspace(0, 100, 101))
    return model


def create_trichloroethylene(parameter_values=None):
    """
    UNCA iGEM 2017 Metabolic Channel.
    """

    # initialize Model
    model = Model(name="Trichloroethylene")

    # Species
    A = Species(name='TCE', initial_value=300)
    B = Species(name='Epoxide', initial_value=120)
    C = Species(name='Dichloracatate', initial_value=0)
    D = Species(name='LossOfOneCL', initial_value=0)
    E = Species(name='Glyoxylate', initial_value=0)
    F = Species(name='Output', initial_value=0)
    model.add_species([A, B, C, D, E, F])

    # Parameters
    K1 = Parameter(name='K1', expression=0.00045 * 0.000025)
    K2 = Parameter(name='K2', expression=14)
    K3 = Parameter(name='K3', expression=0.033)
    K4 = Parameter(name='K4', expression=500 * 0.0001)
    K5 = Parameter(name='K5', expression=500 * 0.0001)
    model.add_parameter([K1, K2, K3, K4, K5])

    # Reactions
    J1 = Reaction(name="J1", reactants={A: 1}, products={B: 1}, rate=K2)

    J2 = Reaction(name="J2", reactants={B: 1}, products={C: 1}, rate=K3)

    J3 = Reaction(name="J3", reactants={C: 1}, products={D: 1}, rate=K4)

    J4 = Reaction(name="J4", reactants={D: 1}, products={E: 1}, rate=K4)

    J5 = Reaction(name="J5", reactants={E: 1}, products={F: 1}, rate=K5)
    model.add_reaction([J1, J2, J3, J4, J5])
    model.timespan(np.linspace(0, 10000, 100))
    return model


def create_lac_operon(parameter_values=None):
    """
    Heath LS, Cao Y. Problem Solving Handbook in Computational Biology and Bioinformatics. Springer; 2014.
    """

    # initialize Model
    model = Model(name="LacOperon")
    
    # Species
    s1 = Species(name='MR', initial_value=0)
    s2 = Species(name='R', initial_value=0)
    s3 = Species(name='R2', initial_value=0)
    s4 = Species(name='O', initial_value=1)
    s5 = Species(name='R2O', initial_value=0)
    s6 = Species(name='I', initial_value=0)
    s7 = Species(name='Iex', initial_value=0)
    s8 = Species(name='I2R2', initial_value=0)
    s9 = Species(name='MY', initial_value=0)
    s10 = Species(name='Y', initial_value=0)
    s11 = Species(name='YIex', initial_value=0)
    s12 = Species(name='Ytot', initial_value=0)
    
    model.add_species(
        [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12])

    # Parameters
    k1 = Parameter(name='k1', expression=0.111)
    k2 = Parameter(name='k2', expression=15.0)
    k3 = Parameter(name='k3', expression=103.8)
    k4 = Parameter(name='k4', expression=0.001)
    k5 = Parameter(name='k5', expression=1992.7)
    k6 = Parameter(name='k6', expression=2.40)
    k7 = Parameter(name='k7', expression=1.293e-6)
    k8 = Parameter(name='k8', expression=12.0)
    k9 = Parameter(name='k9', expression=1.293e-6)
    k10 = Parameter(name='k10', expression=9963.2)
    k11 = Parameter(name='k11', expression=0.50)
    k12 = Parameter(name='k12', expression=0.010)
    k13 = Parameter(name='k13', expression=30)
    k14 = Parameter(name='k14', expression=0.249)
    k15 = Parameter(name='k15', expression=0.10)
    k16 = Parameter(name='k16', expression=60000)
    k17 = Parameter(name='k17', expression=0.920)
    k18 = Parameter(name='k18', expression=0.920)
    k19 = Parameter(name='k19', expression=0.462)
    k20 = Parameter(name='k20', expression=0.462)
    k21 = Parameter(name='k21', expression=0.20)
    k22 = Parameter(name='k22', expression=0.20)
    k23 = Parameter(name='k23', expression=0.20)
    k24 = Parameter(name='k24', expression=0.20)
    k25 = Parameter(name='k25', expression=0.20)
    model.add_parameter(
        [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25])

    # Reactions
    j1 = Reaction(name="j1", reactants={s12: 1}, products={s1: 1}, rate=k1)
    j2 = Reaction(name="j2", reactants={s1: 1}, products={s1: 1, s2: 1}, rate=k2)
    j3 = Reaction(name="j3", reactants={s2: 2}, products={s3: 1}, rate=k3)
    j4 = Reaction(name="j4", reactants={s3: 1}, products={s2: 2}, rate=k4)
    j5 = Reaction(name="j5", reactants={s3: 1, s4: 1}, products={s5: 1}, rate=k5)
    j6 = Reaction(name="j6", reactants={s5: 1}, products={s3: 1, s4: 1}, rate=k6)
    j7 = Reaction(name="j7", reactants={s6: 2, s3: 1}, products={s8: 1}, propensity_function='((3e-7)/((8e-16)*(6.0221367e14))**(2))*(R2)*(I)*(I-1)')
    j8 = Reaction(name="j8", reactants={s8: 1}, products={s6: 2, s3: 1}, rate=k8)
    j9 = Reaction(name="j9", reactants={s6: 2, s5: 1}, products={s8: 1, s4: 1}, propensity_function='((3e-7)/((8e-16)*(6.0221367e14))**(2))*(R2O)*(I)*(I-1)')
    j10 = Reaction(name="j10", reactants={s8: 1, s4: 1}, products={s6: 2, s5: 1}, rate=k10)
    j11 = Reaction(name="j11", reactants={s4: 1}, products={s4: 1, s9: 1}, rate=k11)
    j12 = Reaction(name="j12", reactants={s5: 1}, products={s5: 1, s9: 1}, rate=k12)
    j13 = Reaction(name="j13", reactants={s9: 1}, products={s9: 1, s10: 1}, rate=k13)
    j14 = Reaction(name="j14", reactants={s10: 1, s7: 1}, products={s11: 1}, rate=k14)
    j15 = Reaction(name="j15", reactants={s11: 1}, products={s10: 1, s7: 1}, rate=k15)
    j16 = Reaction(name="j16", reactants={s11: 1}, products={s10: 1, s6: 1}, rate=k16)
    j17 = Reaction(name="j17", reactants={s7: 1}, products={s6: 1}, rate=k17)
    j18 = Reaction(name="j18", reactants={s6: 1}, products={s7: 1}, rate=k18)
    j19 = Reaction(name="j19", reactants={s1: 1}, products={s12: 1}, rate=k19)
    j20 = Reaction(name="j20", reactants={s9: 1}, products={s12: 1}, rate=k20)
    j21 = Reaction(name="j21", reactants={s2: 1}, products={s12: 1}, rate=k21)
    j22 = Reaction(name="j22", reactants={s3: 1}, products={s12: 1}, rate=k22)
    j23 = Reaction(name="j23", reactants={s10: 1}, products={s12: 1}, rate=k23)
    j24 = Reaction(name="j24", reactants={s11: 1}, products={s6: 1}, rate=k24)
    j25 = Reaction(name="j25", reactants={s8: 1}, products={s6: 2}, rate=k25)
    model.add_reaction(
        [j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17, j18, j19, j20, j21, j22, j23, j24, j25])
    model.timespan(np.linspace(0, 100, 10))
    return model


def create_schlogl(parameter_values=None):
    """
    Schlogl F. Chemical reaction models for non-equilibrium phase transitions. Zeitschrift for Physik.
    1972;253: 147–161. doi:10.1007/bf01379769
    """

    # initialize Model
    model = Model(name="Schlogl")

    # Species
    s1 = Species(name='A', initial_value=300)
    s2 = Species(name='B', initial_value=300)
    s3 = Species(name='C', initial_value=300)
    s4 = Species(name='X', initial_value=300)

    model.add_species([s1, s2, s3, s4])

    k1 = Parameter(name='k1', expression=1)
    k2 = Parameter(name='k2', expression=1)

    model.add_parameter([k1, k2])

    j1 = Reaction(name="j1", reactants={s1: 1, s4: 1}, products={s4: 2}, rate=k1)
    j2 = Reaction(name="j2", reactants={s2: 1, s4: 1}, products={s3: 1}, rate=k2)

    model.add_reaction([j1, j2])
    model.timespan(np.linspace(0, 100000, 100))
    return model


def create_michaelis_menten(parameter_values=None):
    # initialize Model
    model = Model(name="Michaelis_Menten")

    # parameters
    rate1 = Parameter(name='rate1', expression=0.0017)
    rate2 = Parameter(name='rate2', expression=0.5)
    rate3 = Parameter(name='rate3', expression=0.1)
    model.add_parameter([rate1, rate2, rate3])

    # Species
    A = Species(name='A', initial_value=301)
    B = Species(name='B', initial_value=120)
    C = Species(name='C', initial_value=0)
    D = Species(name='D', initial_value=0)
    model.add_species([A, B, C, D])

    # reactions
    r1 = Reaction(name="r1", reactants={A: 1, B: 1}, products={C: 1}, rate=rate1)

    r2 = Reaction(name="r2", reactants={C: 1}, products={A: 1, B: 1}, rate=rate2)

    r3 = Reaction(name="r3", reactants={C: 1}, products={B: 1, D: 1}, rate=rate3)
    model.add_reaction([r1, r2, r3])
    model.timespan(np.linspace(0, 100, 101))
    return model


def create_toggle_switch(parameter_values=None):
    """
    Gardner et al. Nature (1999)Construction of a genetic toggle switch in Escherichia coli
    (Transcription from 
    """

    # Initialize the model.
    model = Model(name="Toggle_Switch")
    # Species
    A = Species(name='A', initial_value=10)
    B = Species(name='B', initial_value=10)
    model.add_species([A, B])
    # Parameters
    alpha1 = Parameter(name='alpha1', expression=10)
    alpha2 = Parameter(name='alpha2', expression=10)
    beta = Parameter(name='beta', expression=2)
    gamma = Parameter(name='gamma', expression=2)
    mu = Parameter(name='mu', expression=1)
    model.add_parameter([alpha1, alpha2, beta, gamma, mu])
    # Reactions
    model.add_reaction(Reaction(name="cu", reactants={}, products={'A': 1}, propensity_function="alpha1/(1+pow(B, beta))"))
    model.add_reaction(Reaction(name="cv", reactants={}, products={'B': 1}, propensity_function="alpha2/(1+pow(A, gamma))"))
    model.add_reaction(Reaction(name="du", reactants={'A': 1}, products={}, rate=model.listOfParameters["mu"]))
    model.add_reaction(Reaction(name="dv", reactants={'B': 1}, products={}, rate=model.listOfParameters["mu"]))
    model.timespan(np.linspace(0, 250, 251))
    return model


def create_decay(parameter_values=None):
    """
    This is a simple example for mass-action degradation of species S.
    """

    # Initialize the model.
    model = Model(name="Example")
    # Species
    S = Species(name='Sp', initial_value=100)
    model.add_species([S])
    # Parameters
    k1 = Parameter(name='k1', expression=3.0)
    model.add_parameter([k1])
    # Reactions
    rxn1 = Reaction(name='S degradation', reactants={S: 1}, products={}, rate=k1)
    model.add_reaction([rxn1])
    model.timespan(np.linspace(0, 20, 101))
    return model


def create_decay_no_tspan(parameter_values=None):
    """
    This is a simple example for mass-action degradation of species S.
    """

    # Initialize the model.
    model = Model(name="Example")
    # Species
    S = Species(name='Sp', initial_value=100)
    model.add_species([S])
    # Parameters
    k1 = Parameter(name='k1', expression=3.0)
    model.add_parameter([k1])
    # Reactions
    rxn1 = Reaction(name='S degradation', reactants={S: 1}, products={}, rate=k1)
    model.add_reaction([rxn1])
    return model


def create_tyson_2_state_oscillator(parameter_values=None):
    """
    Here, as a test case, we run a simple two-state oscillator (Novak & Tyson
    2008) as an example of a stochastic reaction system.
    """

    model = Model(name="tyson-2-state", volume=300.0)

    # Species
    X = Species(name='X', initial_value=int(0.65609071 * 300.0))
    Y = Species(name='Y', initial_value=int(0.85088331 * 300.0))
    model.add_species([X, Y])

    P = Parameter(name='p', expression=2.0)
    kt = Parameter(name='kt', expression=20.0)
    kd = Parameter(name='kd', expression=1.0)
    a0 = Parameter(name='a0', expression=0.005)
    a1 = Parameter(name='a1', expression=0.05)
    a2 = Parameter(name='a2', expression=0.1)
    kdx = Parameter(name='kdx', expression=1.0)
    model.add_parameter([P, kt, kd, a0, a1, a2, kdx])

    # creation of X:
    rxn1 = Reaction(name='X production', reactants={}, products={X: 1},
                    propensity_function='300 * 1.0 / (1.0 + (Y * Y / (300 * 300)))')

    # degradadation of X:
    rxn2 = Reaction(name='X degradation', reactants={X: 1}, products={}, rate=kdx)

    # creation of Y:
    rxn3 = Reaction(name='Y production', reactants={X: 1}, products={X: 1, Y: 1}, rate=kt)

    # degradation of Y:
    rxn4 = Reaction(name='Y degradation', reactants={Y: 1}, products={}, rate=kd)

    # nonlinear Y term:
    rxn5 = Reaction(name='Y nonlin', reactants={Y: 1}, products={}, 
            propensity_function='Y / a0 + a1 * (Y / 300) + a2 * Y * Y / (300 * 300)')

    model.add_reaction([rxn1, rxn2, rxn3, rxn4, rxn5])
    model.timespan(np.linspace(0, 100, 101))
    return model

def create_oregonator(parameter_values = None):
    # Oregonator system
    # http://www.scholarpedia.org/article/Oregonator

    # Superclass initialization
    model = Model(name="Oregonator")

    # Species
    F = Species(name = "F", initial_value = 2)
    A = Species(name = "A", initial_value = 250)
    B = Species(name = "B", initial_value = 500)
    C = Species(name = "C", initial_value = 1000)
    P = Species(name = "P", initial_value = 0)
    model.add_species([F, A, B, C, P])

    # Parameters (rates)
    k1 = Parameter(name = "k1", expression = 2.0)
    k2 = Parameter(name = "k2", expression = 0.1)
    k3 = Parameter(name = "k3", expression = 104)
    k4 = Parameter(name = "k4", expression = 4e-7)
    k5 = Parameter(name = "k5", expression = 26.0)
    model.add_parameter([k1, k2, k3, k4, k5])
      
    # Reactions
    reaction1 = Reaction(name = "reaction1", 
                                    reactants = {B: 1, F: 1}, 
                                    products = {A: 1, F: 1}, 
                                    rate = k1)
    reaction2 = Reaction(name = "reaction2", 
                                    reactants = {A: 1, B: 1}, 
                                    products = {P: 1}, 
                                    rate = k2)
    reaction3 = Reaction(name = "reaction3", 
                                    reactants = {A: 1, F: 1}, 
                                    products = {A: 2, C: 1, F: 1}, 
                                    rate = k3)
    reaction4 = Reaction(name = "reaction4", 
                                    reactants = {A: 2}, 
                                    products = {P: 1}, 
                                    rate = k4)
    reaction5 = Reaction(name = "reaction5", 
                                    reactants = {C: 1, F: 1}, 
                                    products = {B: 1, F: 1}, 
                                    rate = k5)
    model.add_reaction([reaction1, reaction2, reaction3, reaction4, reaction5])
      
    # Set timespan of model
    model.timespan(np.linspace(0, 5, 501))
    return model


def create_degradation(model_name="Degradation", parameter_values=None):
    model = Model(name=model_name)
    
    A = Species(name="A", initial_value=2000)
    model.add_species(A)
    
    k1 = Parameter(name="k1", expression="0.1")
    model.add_parameter(k1)
    
    r1 = Reaction(name="r1", rate="k1", reactants={"A":1})
    model.add_reaction(r1)
    model.timespan(np.arange(0, 151))
    return model


def create_robust_model(parameter_values=None):
    model = Model(name="test1")
    model.volume = 1

    # Parameters
    model.add_parameter(Parameter(name="k1", expression="0.5"))
    model.add_parameter(Parameter(name="k2", expression="3.2e-15"))

    # Variables
    model.add_species(Species(name="s1", initial_value=8000, mode="continuous"))
    model.add_species(Species(name="s2", initial_value=0, mode="continuous"))

    # Reactions
    model.add_reaction(Reaction(name="r1", reactants={'s1': 1}, products={'s2': 1}, rate=model.listOfParameters["k2"]))
    model.add_reaction(Reaction(name="r2", reactants={'s1': 1}, products={}, rate=model.listOfParameters["k1"]))
    model.add_reaction(Reaction(name="r3", reactants={'s1': 2}, products={'s1': 3}, propensity_function="sin(1)"))

    # Event Triggers
    e1_trig = EventTrigger(expression="t > 50", initial_value=True, persistent=True)

    # Event Assignments
    e1_assign_1 = EventAssignment(variable="s1", expression="2000")

    # Events
    model.add_event(Event(name="e1", trigger=e1_trig, assignments=[e1_assign_1], delay="t - 40", priority="0", use_values_from_trigger_time=True))

    # Rate Rules
    model.add_rate_rule(RateRule(name="rr1", formula="sin(0.5)", variable="s1"))

    # Assignment Rules
    model.add_assignment_rule(AssignmentRule(name="rr2", formula="8000/(s1+1)", variable="s2"))

    # Function Definitions
    model.add_function_definition(FunctionDefinition(name="multiply", function="x * y", args=['x', 'y']))

    # Timespan
    model.timespan(np.arange(0, 20, 0.05))
    return model

def create_multi_firing_event(parameter_values=None):                                                                                                                                                                                       
    """                                                                                                                                                                                                                              
    This is a simple example for mass-action degradation of species S.  We will add two events
    to demonstrate the usage of events.  The first event will assign the value '0' to our species
    once time passes 20, and the second event will be triggered once time crosses 30, assigning
    a value of "100" to our species and changing the value of our degradation rate parameter
    "k1" from .01 to .1, causing the species to decay more quickly.
    """                                                                                                                                                                                                                              
                                                                                                                                                                                                                                     
    # Initialize the model.                                                                                                                                                                                                      
    model = Model(name="Example")

    # Species                                                                                                                                                                                                                    
    S = Species(name='Sp', initial_value=100, mode='discrete')
    model.add_species([S])
    
    # Parameters                                                                                                                                                                                                                 
    k1 = Parameter(name='k1', expression=0.01)
    ev_time1 = Parameter(name='ev_time1', expression=20)
    ev_time2 = Parameter(name='ev_time2', expression=30)
    model.add_parameter([k1, ev_time1, ev_time2])    
    
    # Events
    et = EventTrigger(expression='t>ev_time1')
    ea = EventAssignment(variable=S, expression='0')
    ea1 = EventAssignment(variable=ev_time1, expression='ev_time1+15')
    e = Event(name='event1', trigger=et, assignments=[ea, ea1])
    
    et2 = EventTrigger(expression='t>ev_time2', persistent=True)
    ea2 = EventAssignment(variable=S, expression='100')
    ea3 = EventAssignment(variable=k1, expression='.1')
    ea4 = EventAssignment(variable=ev_time2, expression='ev_time2+15')
    e2 = Event(name='event2', trigger=et2, assignments=[ea2, ea3, ea4])
    model.add_event([e, e2])
    
    #Reactions
    r = Reaction(name='R', reactants={S:1}, products={}, rate=k1) #Multiple reactions
    model.add_reaction([r])

    model.timespan(np.linspace(0, 60, 181))
    return model


def create_opioid():
    # Initialize Model
    model = gillespy2.Model(name="Opioid")

    # Define Variables (GillesPy2.Species)
    S = gillespy2.Species(name='Susceptibles', initial_value=1000)
    P = gillespy2.Species(name='Prescribed_Users', initial_value=0)
    A = gillespy2.Species(name='Addicted', initial_value=0)
    R = gillespy2.Species(name='Rehab', initial_value=0)
    Natural_Deaths = gillespy2.Species(name='Natural_Deaths', initial_value=0)
    Addiction_Deaths = gillespy2.Species(name='Addiction_Deaths', initial_value=0)

    # Add Variables to Model
    model.add_species([S,P,A,R,Natural_Deaths,Addiction_Deaths])

    # Define Parameters
    alpha = gillespy2.Parameter(name='alpha', expression= 0.15)
    epsilon = gillespy2.Parameter(name='epsilon', expression= 0.8)
    beta_p = gillespy2.Parameter(name='beta_p', expression= 0.00266)
    beta_a = gillespy2.Parameter(name='beta_a', expression= 0.00094)
    gamma = gillespy2.Parameter(name='gamma', expression= 0.00744)
    zeta = gillespy2.Parameter(name='zeta', expression= 0.2)
    delta = gillespy2.Parameter(name='delta', expression= 0.1)
    sigma = gillespy2.Parameter(name='sigma', expression= 0.9)
    mu = gillespy2.Parameter(name='mu', expression= 0.00729)
    mu_prime = gillespy2.Parameter(name='mu_prime', expression= 0.01159)

    # Add Parameters to Model
    model.add_parameter([alpha, epsilon, beta_p, beta_a, gamma, zeta, delta, sigma, mu, mu_prime])

    # Define Reactions
    SP = gillespy2.Reaction(
        name="SP", reactants={'Susceptibles': 1}, products={'Prescribed_Users': 1}, rate='alpha'
    )
    SA_a = gillespy2.Reaction(name="SA_a", reactants={'Susceptibles': 1}, products={'Addicted': 1}, rate='beta_a')
    SA_p = gillespy2.Reaction(name="SA_p", reactants={'Susceptibles': 1}, products={'Addicted': 1}, rate='beta_p')
    mu_S = gillespy2.Reaction(
        name="mu_S", reactants={'Susceptibles': 1}, products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
    )
    PA = gillespy2.Reaction(name="PA", reactants={'Prescribed_Users': 1}, products={'Addicted': 1}, rate='gamma')
    PS = gillespy2.Reaction(
        name="PS", reactants={'Prescribed_Users': 1}, products={'Susceptibles': 1}, rate='epsilon'
    )
    AR = gillespy2.Reaction(name="AR", reactants={'Addicted': 1}, products={'Rehab': 1}, rate='zeta')
    RA = gillespy2.Reaction(name="RA", reactants={'Rehab': 1}, products={'Addicted': 1}, rate='delta')
    RS = gillespy2.Reaction(name="RS", reactants={'Rehab': 1}, products={'Susceptibles': 1}, rate='sigma')
    mu_P = gillespy2.Reaction(
        name="mu_P", reactants={'Prescribed_Users': 1},
        products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
    )
    mu_R = gillespy2.Reaction(
        name="m_R", reactants={'Rehab': 1}, products={'Susceptibles': 1, 'Natural_Deaths': 1}, rate='mu'
    )
    mu_prime_A = gillespy2.Reaction(
        name="mu_prime_A", reactants={'Addicted': 1},
        products={'Susceptibles': 1, 'Addiction_Deaths': 1}, rate='mu_prime'
    )

    # Add Reactions to Model
    model.add_reaction([SP, PS, SA_a, SA_p, PA, AR, RA, RS, mu_P, mu_R, mu_prime_A, mu_S])

    # Define Timespan
    tspan = gillespy2.TimeSpan.linspace(t=10, num_points=11)

    # Set Model Timespan
    model.timespan(tspan)
    return model


def create_telegraph_model():
    """Creates a stochastic model for the dichotomous Markov process OFF <-> ON"""

    # Intialize the Model with a name of your choosing.
    model = gillespy2.Model(name="on_off")

    """
    Variables (GillesPy2.Species) can be anything that participates in or is produced by a reaction channel.

    - name: A user defined name for the species.
    - initial_value: A value/population count of species at start of simulation.
    """
    ON = gillespy2.Species(name="ON", initial_value=0)
    OFF = gillespy2.Species(name="OFF", initial_value=1)

    # Add the Variables to the Model.
    model.add_species([ON, OFF])

    """
    Parameters are constant values relevant to the system, such as reaction kinetic rates.

    - name: A user defined name for reference.
    - expression: Some constant value.
    """
    kon  = gillespy2.Parameter(name="kon",  expression=0.1)
    koff = gillespy2.Parameter(name="koff", expression=1.0)

    # Add the Parameters to the Model.
    model.add_parameter([kon, koff])

    """
    Reactions are the reaction channels which cause the system to change over time.

    - name: A user defined name for the reaction.
    - reactants: A dictionary with participant reactants as keys, and consumed per reaction as value.
    - products: A dictionary with reaction products as keys, and number formed per reaction as value.
    - rate: A parameter rate constant to be applied to the propensity of this reaction firing.
    - propensity_function: Can be used instead of rate in order to declare a custom propensity function in string format.
    """
    r1 = gillespy2.Reaction(
            name="on",
            reactants={'OFF': 1}, 
            products= {'ON': 1},
            rate='kon'
        )

    r2 = gillespy2.Reaction(
            name="off",
            reactants={'ON': 1}, 
            products= {'OFF': 1},
            rate='koff'
        )


    # Add the Reactions to the Model.
    model.add_reaction([r1, r2])

    # Define the timespan of the model.
    tspan = gillespy2.TimeSpan.linspace(t=100, num_points=100)
    
    # Set the timespan of the Model.
    model.timespan(tspan)
    return model



__all__ = [
    'create_trichloroethylene', 'create_lac_operon', 'create_schlogl',
    'create_michaelis_menten', 'create_decay', 'create_decay_no_tspan',
    'create_tyson_2_state_oscillator', 'create_oregonator', 'create_vilar_oscillator',
    'create_dimerization', 'create_degradation', 'create_robust_model',
    'create_multi_firing_event', 'create_toggle_switch','create_telegraph_model',
]
