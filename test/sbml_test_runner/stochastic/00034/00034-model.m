(*

category:        Test
synopsis:        Basic two-reaction system.
componentTags:   Compartment, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

This model is the same as model 00030, with the 'monomer' species ('P') eliminated.

The model contains:
* 1 species (P2)
* 2 parameters (k1, k2)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Dimerisation: -> P2 | $0.5 * k1 * (100 - 2 * P2) * (99 - 2 * P2)$ |
| Disassociation: P2 -> | $k2 * P2$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species P2 | $0$ | variable |
| Initial value of parameter k1 | $0.001$ | constant |
| Initial value of parameter k2 | $0.01$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
