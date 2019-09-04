@model:3.1.1=BatchImmigrationDeath02 "Batch Immigration-Death (004), variant 02"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=0 s
@parameters
 Alpha=1
 Mu=0.4
@reactions
@r=Immigration
 -> 10X
 Alpha
@r=Death
 X -> 
 Mu*X
