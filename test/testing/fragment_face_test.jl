using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

include("/Users/paoluzzi/.julia/v0.6/LinearAlgebraicRepresentation/test/testing/fragment_face.jl")

V,FV = Lar.sphere()([3,4])
EV = Lar.simplexFacets(FV)
obj = V,FV,EV
V,FV,EV = Lar.struct2lar(Lar.Struct([ obj, Lar.t(0.3,0,0.4), obj ]))
Plasm.view(V,FV)



# >>>>> input >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


V, EV, FE, space_idx = preprocessing(V,EV,FV)
for sigma=1:length(FV)
	nV, nEV, nFE = computefragments(V, EV, FE, space_idx, sigma)
end
#sigma=8

