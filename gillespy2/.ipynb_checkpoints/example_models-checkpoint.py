from gillespy2.gillespy2 import *


class Trichloroethylene(Model):
    """
    UNCA iGEM 2017 Metabolic Channel.
    """

    def __init__(self, parameter_values=None):
        # initialize Model
        Model.__init__(self, name="Trichloroethylene")

        # Species
        A = Species(name='TCE', initial_value=300)
        B = Species(name='Epoxide', initial_value=120)
        C = Species(name='Dichloracatate', initial_value=0)
        D = Species(name='LossOfOneCL', initial_value=0)
        E = Species(name='Glyoxylate', initial_value=0)
        F = Species(name='Output', initial_value=0)
        self.add_species([A, B, C, D, E, F])

        # Parameters
        K1 = Parameter(name='K1', expression=0.00045 * 0.000025)
        K2 = Parameter(name='K2', expression=14)
        K3 = Parameter(name='K3', expression=0.033)
        K4 = Parameter(name='K4', expression=500 * 0.0001)
        K5 = Parameter(name='K5', expression=500 * 0.0001)
        self.add_parameter([K1, K2, K3, K4, K5])

        # Reactions
        J1 = Reaction(name="J1", reactants={A: 1}, products={B: 1}, rate=K2)

        J2 = Reaction(name="J2", reactants={B: 1}, products={C: 1}, rate=K3)

        J3 = Reaction(name="J3", reactants={C: 1}, products={D: 1}, rate=K4)

        J4 = Reaction(name="J4", reactants={D: 1}, products={E: 1}, rate=K4)

        J5 = Reaction(name="J5", reactants={E: 1}, products={F: 1}, rate=K5)
        self.add_reaction([J1, J2, J3, J4, J5])
        self.timespan(numpy.linspace(0, 10000, 100))


class LacOperon(Model):
    """
    Heath LS, Cao Y. Problem Solving Handbook in Computational Biology and Bioinformatics. Springer; 2014.
    """

    def __init__(self, parameter_values=None):
        # initialize Model
        Model.__init__(self, name="LacOperon")

        # Species
        s1 = Species(name='PLac', initial_value=300)
        s2 = Species(name='RNAP', initial_value=120)
        s3 = Species(name='PLacRNAP', initial_value=0)
        s4 = Species(name='TrLacZ1', initial_value=0)
        s5 = Species(name='TrLacZ2', initial_value=0)
        s6 = Species(name='TrLacY1', initial_value=0)
        s7 = Species(name='TrLacY2', initial_value=0)
        s8 = Species(name='RBSLacY', initial_value=300)
        s9 = Species(name='RBSLacZ', initial_value=120)
        s10 = Species(name='Ribosome', initial_value=0)
        s11 = Species(name='RbsRibosomeLacY', initial_value=0)
        s12 = Species(name='RbsRibosomeLacZ', initial_value=0)
        s13 = Species(name='TrRbsLacZ', initial_value=0)
        s14 = Species(name='TrRbsLacY', initial_value=0)
        s15 = Species(name='LacZ', initial_value=0)
        s16 = Species(name='LacY', initial_value=0)
        s17 = Species(name='dgrLacZ', initial_value=0)
        s18 = Species(name='dgrLacY', initial_value=0)
        s19 = Species(name='dgrRbsLacZ', initial_value=0)
        s20 = Species(name='dgrRbsLacY', initial_value=0)
        s21 = Species(name='lactose', initial_value=0)
        s22 = Species(name='LacZlactose', initial_value=0)
        s23 = Species(name='product', initial_value=0)

        self.add_species(
            [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22,
             s23])

        # Parameters
        k1 = Parameter(name='k1', expression=0.17)
        k2 = Parameter(name='k2', expression=10)
        k3 = Parameter(name='k3', expression=1)
        k4 = Parameter(name='k4', expression=1)
        k5 = Parameter(name='k5', expression=0.015)
        k6 = Parameter(name='k6', expression=1)
        k7 = Parameter(name='k7', expression=0.36)
        k8 = Parameter(name='k8', expression=0.17)
        k9 = Parameter(name='k9', expression=0.17)
        k10 = Parameter(name='k10', expression=0.45)
        k11 = Parameter(name='k11', expression=0.45)
        k12 = Parameter(name='k12', expression=0.4)
        k13 = Parameter(name='k13', expression=0.4)
        k14 = Parameter(name='k14', expression=0.015)
        k15 = Parameter(name='k15', expression=0.036)
        k16 = Parameter(name='k16', expression=6.42e-5)
        k17 = Parameter(name='k17', expression=6.42e-5)
        k18 = Parameter(name='k18', expression=0.3)
        k19 = Parameter(name='k19', expression=0.3)
        k20 = Parameter(name='k20', expression=9.52e-5)
        k21 = Parameter(name='k21', expression=431)
        k22 = Parameter(name='22', expression=14)
        self.add_parameter(
            [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22])

        # Reactions
        j1 = Reaction(name="j1", reactants={s1: 1, s2: 2}, products={s3: 1}, rate=k1)
        j2 = Reaction(name="j2", reactants={s3: 1}, products={s1: 1, s2: 1}, rate=k2)
        j3 = Reaction(name="j3", reactants={s3: 1}, products={s4: 1}, rate=k3)
        j4 = Reaction(name="j4", reactants={s4: 1}, products={s9: 1, s1: 1, s5: 1}, rate=k4)
        j5 = Reaction(name="j5", reactants={s5: 1}, products={s6: 1}, rate=k5)
        j6 = Reaction(name="j6", reactants={s6: 1}, products={s8: 1, s7: 1}, rate=k6)
        j7 = Reaction(name="j7", reactants={s7: 1}, products={s2: 1}, rate=k7)
        j8 = Reaction(name="j8", reactants={s10: 1, s9: 1}, products={s12: 1}, rate=k8)
        j9 = Reaction(name="j9", reactants={s10: 1, s8: 1}, products={s11: 1}, rate=k9)
        j10 = Reaction(name="j10", reactants={s12: 1}, products={s10: 1, s9: 1}, rate=k10)
        j11 = Reaction(name="j11", reactants={s11: 1}, products={s10: 1, s8: 1}, rate=k11)
        j12 = Reaction(name="j12", reactants={s12: 1}, products={s13: 1, s9: 1}, rate=k12)
        j13 = Reaction(name="j13", reactants={s11: 1}, products={s14: 1, s8: 1}, rate=k13)
        j14 = Reaction(name="j14", reactants={s13: 1}, products={s15: 1}, rate=k14)
        j15 = Reaction(name="j15", reactants={s14: 1}, products={s16: 1}, rate=k15)
        j16 = Reaction(name="j16", reactants={s15: 1}, products={s17: 1}, rate=k16)
        j17 = Reaction(name="j17", reactants={s16: 1}, products={s18: 1}, rate=k17)
        j18 = Reaction(name="j18", reactants={s9: 1}, products={s19: 1}, rate=k18)
        j19 = Reaction(name="j19", reactants={s8: 1}, products={s20: 1}, rate=k19)
        j20 = Reaction(name="j20", reactants={s15: 1, s21: 1}, products={s22: 1}, rate=k20)
        j21 = Reaction(name="j21", reactants={s22: 1}, products={s15: 1, s23: 1}, rate=k21)
        j22 = Reaction(name="j22", reactants={s16: 1}, products={s21: 1, s16: 1}, rate=k22)
        self.add_reaction(
            [j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14, j15, j16, j17, j18, j19, j20, j21, j22])
        self.timespan(numpy.linspace(0, 10000, 100))


class Schlogl(Model):
    """
    Schlogl F. Chemical reaction models for non-equilibrium phase transitions. Zeitschrift for Physik.
    1972;253: 147–161. doi:10.1007/bf01379769
    """

    def __init__(self, parameter_values=None):
        # initialize Model
        Model.__init__(self, name="Schlogl")

        # Species
        s1 = Species(name='A', initial_value=300)
        s2 = Species(name='B', initial_value=300)
        s3 = Species(name='C', initial_value=300)
        s4 = Species(name='X', initial_value=300)

        self.add_species([s1, s2, s3, s4])

        k1 = Parameter(name='k1', expression=1)
        k2 = Parameter(name='k2', expression=1)

        self.add_parameter([k1, k2])

        j1 = Reaction(name="j1", reactants={s1: 1, s4: 1}, products={s4: 2}, rate=k1)
        j2 = Reaction(name="j2", reactants={s2: 1, s4: 1}, products={s3: 1}, rate=k2)

        self.add_reaction([j1, j2])
        self.timespan(numpy.linspace(0, 10000, 100))


class MichaelisMenten(Model):
    def __init__(self, parameter_values=None):
        # initialize Model
        Model.__init__(self, name="Michaelis_Menten")

        # parameters
        rate1 = Parameter(name='rate1', expression=0.0017)
        rate2 = Parameter(name='rate2', expression=0.5)
        rate3 = Parameter(name='rate3', expression=0.1)
        self.add_parameter([rate1, rate2, rate3])

        # Species
        A = Species(name='A', initial_value=301)
        B = Species(name='B', initial_value=120)
        C = Species(name='C', initial_value=0)
        D = Species(name='D', initial_value=0)
        self.add_species([A, B, C, D])

        # reactions
        r1 = Reaction(name="r1", reactants={A: 1, B: 1}, products={C: 1}, rate=rate1)

        r2 = Reaction(name="r2", reactants={C: 1}, products={A: 1, B: 1}, rate=rate2)

        r3 = Reaction(name="r3", reactants={C: 1}, products={B: 1, D: 1}, rate=rate3)
        self.add_reaction([r1, r2, r3])
        self.timespan(numpy.linspace(0, 100, 101))