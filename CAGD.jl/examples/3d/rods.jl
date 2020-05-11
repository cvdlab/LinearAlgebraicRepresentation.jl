todisplay = VERSION <= VersionNumber("1.2") ? true : false

using QHull
using Random
using SparseArrays
using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
if todisplay
	using ViewerGL
	GL = ViewerGL
	include("../examples/3d/views.jl")
end

V,CV = Lar.rod()([3,1])
ch = QHull.chull(convert(Lar.Points,V'))
FV = ch.simplices
EV = Lar.simplexFacets(FV)
GL.VIEW([ GL.GLGrid(V,EV,GL.COLORS[1],1) ]);

cyl = Lar.Struct([ Lar.t(0,0,0), (V,FV,EV) ]) #(0,0,-1.5)
cyl = Lar.struct2lar(cyl)

#V,FV = ([1.0 1.0 -0.5 -0.5 -0.5 -0.5; 0.0 0.0 0.86603 0.86603 -0.86603 -0.86603; 0.0 3.0 0.0 3.0 0.0 3.0],
#Array{Int64,1}[[3, 4, 1, 2], [3, 4, 5, 6], [5, 6, 1, 2], [3, 1, 5], [2, 4, 6]])
#EV = [[1,3],[3,5],[5,1], [2,4],[4,6],[6,2], [1,2],[3,4],[5,6]]
#cyl = (V,FV,EV)

tris = Lar.Struct([
    cyl,
    Lar.Struct([Lar.r(pi/2,0,0), cyl ]) ,
    Lar.Struct([Lar.r(0,pi/2,0), cyl ]) ,
    Lar.Struct([Lar.r(0,pi,0), cyl ]) ,
])
V,FV,EV = Lar.struct2lar(tris)

model = CAGD.Model(V)
scs = convert(Lar.ChainOp, Lar.coboundary_1(model.G, FV, EV))
CAGD.addModelCells!(model, 1, EV, signed = true)
CAGD.addModelCells!(model, 2, scs)

if todisplay  displayModel(model)  end

split_model = CAGD.pairwise_decomposition(model)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)

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

