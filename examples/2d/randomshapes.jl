using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function randomshapes()
	V,EV = Lar.randomcuboids(100, .4)
	V = GL.normalize2(V,flag=true)
	model2d = V,EV
	Sigma =  Lar.spaceindex(model2d);
	GL.VIEW([ GL.GLLines(V,EV, GL.COLORS[1]) ]);

@show(V,EV)
	W,EW = Lar.fraglines(1.5,1.5,1.5)((V,EV))
@show(W,EW)
	
	GL.VIEW([ GL.GLLines(W,EW, GL.COLORS[1]) ]);

	W = convert(Lar.Points, V')
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)
	V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

	triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	W = convert(Lar.Points, V')
	GL.VIEW(GL.GLExplode(W,FVs,1.2,1.2,1.2,99,1));
	GL.VIEW(GL.GLExplode(W,FVs,1,1,1,99,1));
end

randomshapes()
