using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

V,EV = Lar.randomcuboids(250, .2)
V = Plasm.normalize(V,flag=true)
Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
V = convert(Lar.Points, V')
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
V = convert(Lar.Points, W')
Plasm.viewcolor(V::Lar.Points, EVs::Array{Lar.Cells})

model = V,EVs
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))



