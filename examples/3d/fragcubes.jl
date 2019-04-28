using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
using SparseArrays,LinearAlgebra

V,(VV,EV,FV,CV) = Lar.cuboid([0.5,0.5,0.5],true,[-0.5,-0.5,-0.5])
mybox = (V,CV,FV,EV)

twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/3,0,0), Lar.r(0,0,pi/6), mybox ])
V,CV,FV,EV = Lar.struct2lar(twocubes)
Plasm.view(Plasm.numbering(0.25)((V,[[[k] for k=1:size(V,2)],EV,FV])))

sp_idx = Lar.spaceindex((V,FV))
sigma = 6
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
FE = [findnz(cop_FE[k,:])[1] for k=1:size(cop_FE,1)]
FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]





verts, edges = Lar.fragface(V, EV, FV, FE, sp_idx, sigma)

Plasm.view(verts, edges)
