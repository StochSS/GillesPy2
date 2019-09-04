@model:3.1.1=BirthDeath05 "Birth-death model (001), variant 05"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=10000 s
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Lambda*X
@r=Death
 X -> 
 Mu*X
