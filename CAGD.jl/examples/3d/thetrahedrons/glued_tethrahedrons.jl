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

V = [
    0.0 +1  0  0  0
    0.0  0  0 +1 -1
    0.0  0 +1  0  0
];

W = [
     0.0  0.0  0.0
    +1.0  0.0  0.0
     0.0  0.0 +1.0
     0.0 +1.0  0.0
     0.0 -1.0  0.0
];

EV = SparseArrays.sparse(Int8.([
    -1 +1  0  0  0
    -1  0 +1  0  0
    -1  0  0 +1  0
    -1  0  0  0 +1
     0 -1 +1  0  0
     0 -1  0 +1  0
     0 -1  0  0 +1
     0  0 -1 +1  0
     0  0 -1  0 +1
]));

FE = SparseArrays.sparse(Int8.([
    +1 -1  0  0 +1   0  0  0  0
    +1  0 -1  0  0  +1  0  0  0
    +1  0  0 -1  0   0 +1  0  0
     0 +1 -1  0  0   0  0 +1  0
     0 +1  0 -1  0   0  0  0 +1
     0  0  0  0 +1  -1  0 +1  0
     0  0  0  0 +1   0 -1  0 +1
]));

CF = SparseArrays.sparse(Int8.([
    -1 +1  0 -1  0 +1  0
    +1  0 -1  0 +1  0 -1
     0 -1 +1 +1 -1 -1 +1
]));

model = CAGD.Model(V)
CAGD.addModelCells!(model, 1, EV)
CAGD.addModelCells!(model, 2, FE)
CAGD.addModelCells!(model, 3, CF)



atol = 1e-6;
if todisplay  displayModel(model)  end

split_model = CAGD.facesplitting(model, atol = atol)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergemodel(split_model, err=atol, signed_merge=true)
#split_model = CAGD.facesplitting(congr_model, atol = atol)
#congr_model = CAGD.mergemodel(split_model, err=atol, signed_merge=true)

if todisplay  displayModel(congr_model)  end

gift_model = deepcopy(congr_model)
FCcycles, bicon_comps = CAGD.tgw(congr_model, 3, atol = atol)
CAGD.addModelCells!(gift_model, 3, convert(Lar.ChainOp, FC'))

if todisplay  viewExplode(gift_model, complete = false)  end

# Since first face is inner, then retrieving the external is needed

externals = [bc[CAGD.getExternalBoundary(congr_model, 3, FCcycles, bc)] for bc in bicon_comps]
internals = [setdiff(bicon_comps[i], externals[i]) for i = 1 : length(externals)]
assembled_model = deepcopy(congr_model)

FC = FCcycles[:, internals...]
JK = ones(Int8, length(externals))
FC = [FCcycles * SparseArrays.sparse(externals, JK, JK, size(FCcycles, 2), 1) FC]
CAGD.addModelCells!(assembled_model, 3, convert(Lar.ChainOp, FC'))

if todisplay  viewExplode(assembled_model)  end
