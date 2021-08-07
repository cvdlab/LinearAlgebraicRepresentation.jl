using LinearAlgebraicRepresentation,ViewerGL
Lar = LinearAlgebraicRepresentation
GL = ViewerGL

function randomarrangement2d(V,EV)
	VV = [[k] for k=1:size(V,2)]
	#GL.VIEW( GL.numbering(0.1)((V,[VV, EV])) );

	W = convert(Lar.Points, V')
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)
	V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

	triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	V = convert(Lar.Points, V')
	#GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99));

	W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
	EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
	V = convert(Lar.Points, W')
	#GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,1));

	return V,EVs,FVs
end

V,EV = Lar.randomcuboids(100, .4)
V = GL.normalize2(V,flag=true)

W,EVs,FVs = randomarrangement2d(V,EV)
GL.VIEW(GL.GLExplode(W,EVs,1.2,1.2,1.2,1,1));
GL.VIEW(GL.GLExplode(W,FVs,1.2,1.2,1.2,99,1));
GL.VIEW(GL.GLExplode(W,FVs,1,1,1,99,1));
