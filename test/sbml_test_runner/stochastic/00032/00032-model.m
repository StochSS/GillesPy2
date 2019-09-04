(*

category:        Test
synopsis:        Basic two-reaction system with event.
componentTags:   CSymbolTime, Compartment, EventNoDelay, Parameter, Reaction, Species
testTags:        Amount, HasOnlySubstanceUnits, NonUnityStoichiometry
testType:        StochasticTimeCourse
levels:          2.1, 2.2, 2.3, 2.4, 3.1
generatedBy:     Analytic
packagesPresent: 

This model is the same as model 00030, with an event that resets the two species levels based on time.

The model contains:
* 2 species (P, P2)
* 2 parameters (k1, k2)
* 1 compartment (Cell)

There are 2 reactions:

[{width:30em,margin: 1em auto}|  *Reaction*  |  *Rate*  |
| Dimerisation: 2P -> P2 | $(k1 * P * (P - 1)) / 2$ |
| Disassociation: P2 -> 2P | $k2 * P2$ |]


There is one event:

[{width:30em,margin: 1em auto}|  *Event*  |  *Trigger*  | *Event Assignments* |
| reset | $t >= 25$ | $P = 100$ |
|  |  | $P2 = 0$ |]

The initial conditions are as follows:

[{width:35em,margin: 1em auto}|       | *Value* | *Constant* |
| Initial amount of species P | $100$ | variable |
| Initial amount of species P2 | $0$ | variable |
| Initial value of parameter k1 | $0.001$ | constant |
| Initial value of parameter k2 | $0.01$ | constant |
| Initial volume of compartment 'Cell' | $unknown$ | constant |]

Note: The test data for this model was generated from an analytical
solution of the system of equations.

*)
