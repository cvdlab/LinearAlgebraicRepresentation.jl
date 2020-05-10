todisplay = VERSION <= VersionNumber("1.2") ? true : false

using CAGD
using SparseArrays
using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
if todisplay
    using ViewerGL
	GL = ViewerGL
	include("views.jl")
end

V = 2 * [
    0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0
    0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0
    0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
]
EV = [
    [1, 2], [1, 3], [2, 4], [3, 4],
    [1, 5], [2, 6], [3, 7], [4, 8],
    [5, 6], [5, 7], [6, 8], [7, 8]
]
FE = [
    [1, 2, 3, 4], [9, 10, 11, 12],
    [1, 5, 6, 9], [2, 5, 7, 10], [3, 6, 8, 11], [4, 7, 8, 12]
]
CF = [
    [1, 2, 3, 4, 5, 6]
]

m1 = CAGD.Model(V)
CAGD.addModelCells!(m1, 1, EV, signed = true)
CAGD.addModelCells!(m1, 2, FE, signed = true)
CAGD.addModelCells!(m1, 3, CF, signed = false)

m2 = deepcopy(m1)
m2.G .+= [1.0; 1.0; 1.0]

model = CAGD.uniteModels(m1, m2)

if todisplay  displayModel(model)  end

split_model = CAGD.facesplitting(model)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergemodel(split_model, signed_merge=true)

if todisplay  displayModel(congr_model)  end

gift_model = deepcopy(congr_model)

FCcycles, bicon_comps = CAGD.tgw(congr_model, 3)

externals = [bc[CAGD.getExternalBoundary(congr_model, 3, FCcycles, bc)] for bc in bicon_comps]
internals = [setdiff(bicon_comps[i], externals[i]) for i = 1 : length(externals)]

FC = FCcycles[:, internals...]
JK = ones(Int8, length(externals))
FC = [FCcycles * SparseArrays.sparse(externals, JK, JK, size(FCcycles, 2), 1) FC]

CAGD.addModelCells!(gift_model, 3, convert(Lar.ChainOp, FC'))

if todisplay  viewExplode(gift_model)  end

# Each face is percurred two times
sum(gift_model.T[3], dims = 1)

boolmatrix, arranged_model = CAGD.bool3(model)

arranged_model == gift_model
