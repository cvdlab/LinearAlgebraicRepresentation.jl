# TODO:  TO FIX !!!

using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using SparseArrays, LinearAlgebra
using ViewerGL
GL = ViewerGL

function randomshapes()
	V,EV,FV = Lar.randomcubes(8, .8)
	model3d = V,EV,FV
	VV = [[k] for k=1:size(V,2)]

	Sigma =  Lar.spaceindex((V,FV))
	GL.VIEW( GL.numbering(.15)((V,[VV, EV, FV])) )
	V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components (H & T)
	GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,12,0.99));

	W,bases,coboundaries = Lar.chaincomplex(V,FV,EV)
	EV, FV, CV = bases
	copEV, copFE, copCF = coboundaries # OK

#
EV = Lar.cop2lar(copEV)
FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]
FV = convert(Lar.Cells, FV)
W = convert(Lar.Points, W)
WW = [[k] for k=1:size(W,2)]

GL.VIEW(GL.numbering(.2)((W,[WW,EV]) ));
#




	W = convert(Lar.Points, W')
	triangulated_faces = Lar.triangulate(W, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, Lar.triangulated_faces)


	Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})
end

randomshapes()
