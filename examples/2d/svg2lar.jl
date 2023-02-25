using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL,SparseArrays
GL = ViewerGL

function myshow(filename)
	V, EV = Lar.svg2lar(filename)
	GL.VIEW([ GL.GLLines(V,EV) ])
	return V, EV
end

# myshow("svg/new.svg")
# myshow("svg/curved.svg")
# myshow("svg/twopaths.svg")
# myshow("svg/paths.svg")
# myshow("svg/boundarytest2.svg")
# myshow("svg/tile.svg")
# myshow("svg/interior.svg")
V,EV = myshow("svg/holes.svg")
# V,EV = myshow("svg/Lar.svg")

# subdivision of input edges
W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

# compute containment graph of components
bicon_comps = Lar.Arrangement.biconnected_components(copEV)

# visualization of component graphs
EW = Lar.cop2lar(copEV);
W = convert(Lar.Points, V');
V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components
VV = [[k] for k=1:size(V,2)];
meshes = [ GL.numbering(.1)((V, (VV,EVs[k])), GL.COLORS[k],1) for k=1:length(EVs) ];
GL.VIEW( union(meshes...) );

# final solid visualizations
FE = [SparseArrays.findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
FV = [collect(Set(Lar.cat([EW[e] for e in FE[f]]))) for f=1:length(FE)]
V = convert(Lar.Points, V')
triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
V = convert(Lar.Points, V')
FVs = convert(Array{Lar.Cells}, triangulated_faces)
GL.VIEW(GL.GLExplode(V,FVs,1,1,1,99,1));  

V,FVs,EVs = Lar.arrange2D(V,EW)
GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,99,1));
GL.VIEW(GL.GLExplode(V,FVs,1,1,1,99,1));  

#EVs = Lar.FV2EVs(copEV, copFE) # polygonal face boundaries
#EVs = convert(Array{Array{Array{Int64,1},1},1}, EVs)
#GL.VIEW(GL.GLExplode(V,EVs,1,1,1,3,1));  # TODO:  fix inner holes
