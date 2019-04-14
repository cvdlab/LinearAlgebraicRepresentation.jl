using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm


V,EV = Lar.randomcuboids(100, .2)
V = Plasm.normalize(V,flag=true)
model2d = V,EV

Sigma =  Lar.spaceindex(model2d);
for k=1:length(Sigma) println(k,Sigma[k]) end

Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV])))
Plasm.viewexploded(V,EV)(1.2,1.2,1.2) 	# no numerical errors

W,EW = Lar.fragmentlines((V,EV)) 
Plasm.viewexploded(W,EW)(1.2,1.2,1.2)	
Plasm.view(Plasm.numbering(.015)((W,[[[k] for k=1:size(W,2)], EW])))

V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)
hpcs = [ Plasm.lar2hpc(V,EVs[i]) for i=1:length(EVs) ]
Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EW::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
W = convert(Lar.Points, V')
Plasm.viewcolor(W::Lar.Points, FVs::Array{Lar.Cells})


