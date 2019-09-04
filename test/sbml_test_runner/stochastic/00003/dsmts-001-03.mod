@model:3.1.1=BirthDeath03 "Birth-death model (001), variant 03"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=100 s
@parameters
 Lambda=1
 Mu=1.1
@reactions
@r=Birth
 X ->  2X
 Lambda*X
@r=Death
 X -> 
 Mu*X
