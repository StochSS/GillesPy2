(*

category:        Test
synopsis:        Basic two-reaction system.
componentTags:   Compartment, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits, LocalParameters
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

This model is the same as model 00020, but with local parameters used for both reaction rates, that both shadow a global parameter of the same name.

The model contains:
* 1 species (X)
* 1 parameter (k)
* 2 local parameters (Immigration.k, Death.k)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Immigration: -> X | $k$ |
| Death: X -> | $k * X$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species X | $0$ | variable |
| Initial value of parameter k | $2$ | constant |
| Initial value of local parameter 'Immigration.k' | $1$ | constant |
| Initial value of local parameter 'Death.k' | $0.1$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
