todisplay = VERSION <= VersionNumber("1.2") ? true : false

using CAGD
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

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
model = CAGD.Model(V);
CAGD.addModelCells!(model, 1, convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)))
CAGD.addModelCells!(model, 2, Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells))

# Spatial Arrangement
split_model = CAGD.facesplitting(model)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)

if todisplay  displayModel(congr_model)  end

gift_model = deepcopy(congr_model)
FCcycles, bicon_comps = CAGD.tgw(congr_model, 3)
CAGD.addModelCells!(gift_model, 3, convert(Lar.ChainOp, FCcycles'))

if todisplay  viewExplode(gift_model)  end

externals = [bc[CAGD.getExternalBoundary(congr_model, 3, FCcycles, bc)] for bc in bicon_comps]
internals = [setdiff(bicon_comps[i], externals[i]) for i = 1 : length(externals)]
assembled_model = deepcopy(congr_model)
## To Fix
FC = FCcycles[:, internals...]
JK = ones(Int8, length(externals))
FC = [FCcycles * SparseArrays.sparse(externals, JK, JK, size(FCcycles, 2), 1) FC]
CAGD.addModelCells!(assembled_model, 3, convert(Lar.ChainOp, FC'))

if todisplay  viewExplode(assembled_model)  end