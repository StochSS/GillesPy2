@model:3.1.1=BirthDeath02 "Birth-death model (001), variant 02"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=100 s
@reactions
@r=Birth
 X ->  2X
 Lambda*X : Lambda=0.1
@r=Death
 X -> 
 Mu*X : Mu=0.11
