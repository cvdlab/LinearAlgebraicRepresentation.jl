using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

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
GL.VIEW([ GL.GLGrid(V,FV, GL.Point4d(1,1,1,0.2)) ]);
m = CAGD.Model(V);
CAGD.addModelCells!(m, 1, convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells)))
CAGD.addModelCells!(m, 2, Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells))


sp_idx = CAGD.spaceIndex(m, 2)

de_models = [
	CAGD.face_decomposition(m, face_idx, sp_idx[face_idx])
	for face_idx = 1 : size(m, 2, 1)
]

mfrag = CAGDmergeMultipleModels(de_models)

triangulated_faces = Lar.triangulate(convert(Lar.Points, mfrag.G'), mfrag.T[1:2])
FVs = convert(Array{Lar.Cells}, triangulated_faces);
GL.VIEW( GL.GLExplode(mfrag.G,FVs,expl,expl,expl,99,1) );

mArranged = CAGDdeepcopy(mfrag)
CAGDaddModelCells!(mArranged, 3, Lar.Arrangement.minimal_3cycles(
		convert(Lar.Points, mfrag.G'), mfrag.T[1], mfrag.T[2])
)




V,CVs,FVs,EVs = Lar.pols2tria(mArranged.G, mArranged.T[1], mArranged.T[2], mArranged.T[3])

GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
meshes = GL.GLExplode(V,CVs[1:end],8,4,6,99,0.5);
push!( meshes, GL.GLFrame)
GL.VIEW( meshes );

#=
cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)

EV = Lar.cop2lar(copEV)
FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]
FV = convert(Lar.Cells, FV)
W = convert(Lar.Points, V')
WW = [[k] for k=1:size(W,2)]

GL.VIEW(GL.numbering(.5)((W,[WW,EV]) ));

triangulated_faces = Lar.triangulate(V, [copEV, copFE])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
V = convert(Lar.Points, V')
GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));

EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
=#
