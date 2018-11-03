#using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
L = LinearAlgebraicRepresentation

#include("/Users/paoluzzi/.julia/v0.6/LinearAlgebraicRepresentation/test/testing/fragment_face.jl")

# err = LinearAlgebraicRepresentation.ERR

#V,FV = Lar.sphere()([3,4])
#V,FV = Lar.sphere()([6,9])
#V,FV = Lar.sphere()([12,18]) 
#V,FV = Lar.sphere()([15,24]) 
#V,FV = Lar.sphere()([20,30]) 
#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40])
#V,FV = Lar.sphere()([30,40]) #KO
#V,FV = Lar.sphere()([21,32]) #KO
#V,FV = Lar.sphere()([22,34]) #KO
#V,FV = Lar.sphere()([30,48]) 

# err = 1e-10

#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40]) #KO 

# err = 1e-12

#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40]) #KO 

# err = 1e-7

#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40]) #KO 

# err = 1e-6

#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40]) #KO 

# err = 1e-4

#V,FV = Lar.sphere()([20,40]) 
#V,FV = Lar.sphere()([25,40]) #KO 

# err = 1e-3

#V,FV = Lar.sphere()([20,40]) #  <=== testare con mio merge_vertices: indicazione su errore
#V,FV = Lar.sphere()([25,40]) #KO 


#assembly = L.Struct(repeat([ L.sphere()(), L.t(0.8,0,0), L.r(0,pi/10,0), 	L.r(0,0,pi/5)],8))

assembly = L.Struct(repeat([ L.sphere()([32,48]), L.t(0,1,0), L.r(0,pi/4,0) ],2))
V,FV = Lar.struct2lar(assembly)
EV = Lar.simplexFacets(FV)
Plasm.view(V,FV)
W,FW,EW = V,FV,EV
V,bases,coboundaries = Lar.chaincomplex(V,FV,EV)


