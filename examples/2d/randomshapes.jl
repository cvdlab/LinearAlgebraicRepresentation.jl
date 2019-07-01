using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

function randomshapes()
	V,EV = Lar.randomcuboids(50, .4)
	V = Plasm.normalize(V,flag=true)
	model2d = V,EV
	Sigma =  Lar.spaceindex(model2d);
	GL.VIEW([ GL.GLLines(V,EV, GL.COLORS[1]) ]);

	W,EW = Lar.fraglines(1.5,1.5,1.5)((V,EV))
	GL.VIEW([ GL.GLLines(W,EW, GL.COLORS[1]) ]);

	# # Too slow ...
	# W,EW = Lar.fragmentlines((V,EV))
	# V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)
	# VV = [[k] for k=1:size(V,2)]
	# meshes = [ GL.numbering(0.01)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
	# GL.VIEW( cat(meshes) );

	W = convert(Lar.Points, V')
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)
	V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

	triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	W = convert(Lar.Points, V')
	GL.VIEW(GL.GLExplode(W,FVs,1.2,1.2,1.2,99));
end

randomshapes()
