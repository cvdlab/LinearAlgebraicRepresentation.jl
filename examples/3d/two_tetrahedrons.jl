using Plasm, SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Revise

V, (VV,EV,FV,CV) = Lar.simplex(3, true)
tetra = V, EV,FV,CV
twotetra = Lar.Struct([ tetra, Lar.t(0.25,0.25,0.25), tetra ])
V,EV,FV,CV = Lar.struct2lar(twotetra)
Plasm.view(V,CV)

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

Plasm.view(Plasm.numbering(0.25)((W,[[[k] for k=1:size(W,2)],EV,FV])))

triangulated_faces = Lar.triangulate(V, [copEV, copFE])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
V = convert(Lar.Points, V')
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
model = V,EVs
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
#
#
#  @test EVs==[[[1, 2], [1, 3], [2, 5], [3, 4], [5, 6], [4, 6]],
#  [[4, 5], [5, 6], [4, 6]]                        ,
#  [[1, 2], [1, 7], [2, 7]]                        ,
#  [[1, 3], [1, 7], [3, 7]]                        ,
#  [[2, 5], [3, 4], [2, 7], [3, 7], [5, 8], [4, 8]],
#  [[4, 5], [5, 8], [4, 8]]                        ,
#  [[4, 5], [5, 9], [4, 10], [9, 10]]              ,
#  [[5, 6], [5, 8], [6, 8]]                        ,
#  [[5, 8], [5, 9], [8, 11], [9, 11]]              ,
#  [[4, 6], [4, 8], [6, 8]]                        ,
#  [[4, 8], [4, 10], [8, 11], [10, 11]]            ,
#  [[9, 10], [9, 11], [10, 11]]]
