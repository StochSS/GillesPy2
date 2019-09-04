@model:3.1.1=ImmigrationDeath09 "Immigration-Death (002), variant 09"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=0 s
@parameters
 Alpha=1
 Mu=0.1
@reactions
@r=Immigration
 -> X
 Alpha
@r=Death
 X -> 
 Mu*X
@events
 reset= t>=25 : X=50
