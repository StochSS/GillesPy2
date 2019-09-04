@model:3.1.1=BirthDeath08 "Birth-death model (001), variant 08"
 s=item,t=second,v=litre
@compartments
 Cell=1
@species
 Cell:X=100 s
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
