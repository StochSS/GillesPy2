@model:3.1.1=BirthDeath06 "Birth-death model (001), variant 06"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:X=100 s
 Cell:Sink=0 sb
@parameters
 Lambda=0.1
 Mu=0.11
@reactions
@r=Birth
 X ->  2X
 Lambda*X
@r=Death
 X -> Sink
 Mu*X
