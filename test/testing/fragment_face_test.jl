using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

#include("/Users/paoluzzi/.julia/v0.6/LinearAlgebraicRepresentation/test/testing/fragment_face.jl")

#V,FV = Lar.sphere()([3,4])
#V,FV = Lar.sphere()([6,9])
V,FV = Lar.sphere()([12,18]) #KO
#V,FV = Lar.sphere()([15,24]) #KO
EV = Lar.simplexFacets(FV)
obj = V,FV,EV
V,FV,EV = Lar.struct2lar(Lar.Struct([ obj, Lar.t(0.3,0,0.4), obj ]))
Plasm.view(V,FV)
W,FW,EW = V,FV,EV
V,bases,coboundaries = Lar.chaincomplex(V,FV,EV)


