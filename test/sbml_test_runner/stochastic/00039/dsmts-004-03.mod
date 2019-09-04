@model:3.1.1=BatchImmigrationDeath03 "Batch Immigration-Death (004), variant 03"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=0 s
@parameters
 Alpha=1
 Mu=4
@reactions
@r=Immigration
 -> 100X
 Alpha
@r=Death
 X -> 
 Mu*X
