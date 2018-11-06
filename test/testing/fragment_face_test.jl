#using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
L = LinearAlgebraicRepresentation

#assembly = L.Struct(repeat([ L.sphere()([16,22]), L.t(0.8,0,0), L.r(0,pi/10,0), 	L.r(0,0,pi/5)],8))

assembly = L.Struct(repeat([ L.sphere()([30,48]), L.t(0,1,0), L.r(0,pi/2,0) ],2))
V,FV = Lar.struct2lar(assembly)
EV = Lar.simplexFacets(FV)
Plasm.view(V,FV)
W,FW,EW = V,FV,EV
V,bases,coboundaries = Lar.chaincomplex(V,FV,EV)


