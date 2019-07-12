using SparseArrays
using ViewerGL
GL = ViewerGL
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation


# input generation
V,EV = Lar.randomcuboids(60, .35)
V = GL.normalize2(V,flag=true)
VV = [[k] for k=1:size(V,2)]
GL.VIEW( GL.numbering(.05)((V,[VV, EV]),GL.COLORS[1]) )

# subdivision of input edges
W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

# compute containment graph of components
bicon_comps = Lar.Arrangement.biconnected_components(copEV)

# visualization of component graphs
EW = Lar.cop2lar(copEV)
W = convert(Lar.Points, V')
comps = [ GL.GLLines(W,EW[comp],GL.COLORS[(k-1)%12+1]) for (k,comp) in enumerate(bicon_comps) ]
GL.VIEW(comps);

# visualization of numbered arrangement
VV = [[k] for k=1:size(V,1)]
EV = Lar.cop2lar(copEV);
FE = Lar.cop2lar(copFE);
FV = cat([[EV[e] for e in face] for face in FE])

# final solid visualization
triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
V = convert(Lar.Points, V')
FVs = convert(Array{Lar.Cells}, triangulated_faces)
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99,1));

# polygonal face boundaries
EVs = Lar.FV2EVs(copEV, copFE)
EVs = convert(Array{Array{Array{Int64,1},1},1}, EVs)
GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,1,1));
