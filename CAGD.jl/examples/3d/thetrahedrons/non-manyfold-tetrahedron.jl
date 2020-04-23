todisplay = VERSION <= VersionNumber("1.2") ? true : false

using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD

if todisplay
	using ViewerGL
	GL = ViewerGL
	include("../examples/3d/views.jl")
end

V1, (VV,EV,FV,CV) = Lar.simplex(3, true)
tetra1 = V1, EV, FV, CV
V2 = [
    0.0 0.5 0.1 0.1
    0.0 0.1 0.5 0.1
    0.0 0.1 0.1 0.5
]
tetra2 = V2, EV, FV, CV
twotetra = Lar.Struct([ tetra1, tetra2 ])
V,EV,FV,CV = Lar.struct2lar(twotetra)
# GL.VIEW([ GL.GLGrid(V,FV, GL.Point4d(1,1,1,0.2)) ]);
model = CAGD.Model(V);
CAGD.addModelCells!(model, 1, convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)))
CAGD.addModelCells!(model, 2, Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells))

# Spatial Arrangement
split_model = CAGD.pairwise_decomposition(model)
congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)

gift_model, bicon_comps = CAGD.tgw(congr_model, 3)
