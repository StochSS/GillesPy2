@model:3.1.1=Dimerisation05 "Dimerisation model (003), variant 05"
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
 0.5*k1*(100-2*P2)*(99-2*P2)
@r=Disassociation
 P2 -> 
 k2*P2
