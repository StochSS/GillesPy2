(*

category:        Test
synopsis:        Basic two-reaction system
componentTags:   Compartment, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits, NonUnityStoichiometry
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

A basic two-reaction system.

The model contains:
* 1 species (X)
* 2 parameters (Lambda, Mu)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Birth: X -> 2X | $Lambda * X$ |
| Death: X -> | $Mu * X$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species X | $100$ | variable |
| Initial value of parameter Lambda | $0.1$ | constant |
| Initial value of parameter Mu | $0.11$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
