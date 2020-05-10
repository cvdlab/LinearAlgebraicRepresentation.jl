todisplay = VERSION <= VersionNumber("1.2") ? true : false

using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
if todisplay
    using ViewerGL
    GL = ViewerGL
	include("../views.jl")
end

V = [
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

m1 = CAGD.Model(V)
CAGD.addModelCells!(m1, 1, EV, signed = true)
CAGD.addModelCells!(m1, 2, FE, signed = true)

rx = Lar.r(π/6,0,0)
ry = Lar.r(0,π/3,0)
rot = rx * ry
m2 = deepcopy(m1)
m2.G = rot[1:3, :] * [m2.G; ones(1, size(m2, 0, 2))]

model = CAGD.uniteModels(m1, m2)

if todisplay  displayModel(model)  end

split_model = CAGD.facesplitting(model)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergemodel(split_model, signed_merge=true)

if todisplay  displayModel(congr_model)  end

FC, bicon_comps = CAGD.tgw(congr_model, 3)
gift_model = deepcopy(congr_model);
CAGD.addModelCells!(gift_model, 3, convert(Lar.ChainOp, FC'));

if todisplay  viewExplode(gift_model)  end
