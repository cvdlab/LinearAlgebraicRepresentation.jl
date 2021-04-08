using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

include("./twocubegrids.jl")

# generation of the chain 3-complex of two intersecting cell 2-complexes in 3-space
W, copEV, copFE, copCF = twocubegrids(2,2,2);
# generation of arrays of simplicial complexes: by body, by face, and by edges
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)
# generation and explosion of meshes by body
meshes = GL.GLExplode(V,CVs[2:end],3,3,3,99,.5)
GL.VIEW(meshes);
meshes = GL.GLExplode(V,FVs[2:end],3,3,3,99,.5)
# visualization
GL.VIEW(meshes);
