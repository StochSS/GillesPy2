@model:3.1.1=Dimerisation02 "Dimerisation model (003), variant 02"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:P=1000 s
 Cell:P2=0 s
@parameters
 k1=0.0002
 k2=0.004
@reactions
@r=Dimerisation
 2P ->  P2
 k1*P*(P-1)/2
@r=Disassociation
 P2 -> 2P 
 k2*P2
