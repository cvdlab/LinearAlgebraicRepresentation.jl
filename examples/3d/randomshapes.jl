using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, SparseArrays, LinearAlgebra

function randomshapes()
	V,EV,FV = Lar.randomcubes(5, .8)
	V = Plasm.normalize3D(V; flag=true)
	model3d = V,EV,FV

	Sigma =  Lar.spaceindex((V,FV));

	Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV])))
	Plasm.viewexploded(V,EV)(1.2,1.2,1.2) 	# no numerical errors

	V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components (H & T)
	hpcs = [ Plasm.lar2hpc(V,EVs[i]) for i=1:length(EVs) ]
	Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

	W,bases,coboundaries = Lar.chaincomplex(V,FV,EV)
	EV, FV, CV = bases
	copEV, copFE, copCF = coboundaries

	W = convert(Lar.Points, W')
	triangulated_faces = Lar.triangulate(W, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})
end

randomshapes()
