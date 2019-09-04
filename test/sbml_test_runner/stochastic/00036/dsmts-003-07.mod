@model:3.1.1=Dimerisation07 "Dimerisation model (003), variant 07"
 s=item,t=second,v=litre
@compartments
 Cell
@species
 Cell:P2=0 s
@parameters
 k1=0.001
 k2=0.01
@reactions
@r=Dimerisation
 ->  P2
 k1*(100-2*P2)*((100-2*P2)-1)/2
@r=Disassociation
 P2 -> 
 k2*P2
