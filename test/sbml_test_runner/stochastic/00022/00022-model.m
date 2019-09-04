(*

category:        Test
synopsis:        Basic two-reaction system.
componentTags:   Compartment, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits, LocalParameters
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

This model is the same as case 00020, but with a different 'Immigration' rate that uses a local parameter that shadows a global parameter of the same name.

The model contains:
* 1 species (X)
* 2 parameters (Alpha, Mu)
* 1 local parameter (Immigration.Alpha)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Immigration: -> X | $Alpha$ |
| Death: X -> | $Mu * X$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species X | $0$ | variable |
| Initial value of parameter Alpha | $10$ | constant |
| Initial value of parameter Mu | $0.1$ | constant |
| Initial value of local parameter 'Immigration.Alpha' | $5$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
