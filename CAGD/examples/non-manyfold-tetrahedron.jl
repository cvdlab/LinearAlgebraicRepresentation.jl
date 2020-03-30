using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL
CAGD = CAGD

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
m = CAGD.Model(V);
CAGD.addModelCells!(m, 1, convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)))
CAGD.addModelCells!(m, 2, Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells))

# Spatial Arrangement

## Pairwise Decomposition

# split_model = CAGD.pairwise_decomposition(m)

sp_idx = CAGD.spaceIndex(m, 2)

de_models = [
	CAGD.face_decomposition(m, face_idx, sp_idx[face_idx])
	for face_idx = 1 : size(m, 2, 1)
]

split_model = CAGD.uniteMultipleModels(de_models)

## Congruence via CCE algorithm

congr_model = CAGD.mergeModelVertices(split_model, signed=true)
