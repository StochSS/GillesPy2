(*

category:        Test
synopsis:        Basic two-reaction system with event.
componentTags:   CSymbolTime, Compartment, EventNoDelay, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

This model is based off of model 00020, but with an event halfway through the time course.

The model contains:
* 1 species (X)
* 2 parameters (Alpha, Mu)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Immigration: -> X | $Alpha$ |
| Death: X -> | $Mu * X$ |]


There is one event:

[{width:30em,margin: 1em auto}|  *Event*  |  *Trigger*  | *Event Assignments* |
| reset | $t >= 25$ | $X = 50$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species X | $0$ | variable |
| Initial value of parameter Alpha | $1$ | constant |
| Initial value of parameter Mu | $0.1$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
