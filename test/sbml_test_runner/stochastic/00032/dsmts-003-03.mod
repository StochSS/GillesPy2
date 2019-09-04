@model:3.1.1=Dimerisation03 "Dimerisation model (003), variant 03"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:P=100 s
 Cell:P2=0 s
@parameters
 k1=0.001
 k2=0.01
@reactions
@r=Dimerisation
 2P ->  P2
 k1*P*(P-1)/2
@r=Disassociation
 P2 -> 2P 
 k2*P2
@events
 reset= t>=25 : P=100; P2=0
