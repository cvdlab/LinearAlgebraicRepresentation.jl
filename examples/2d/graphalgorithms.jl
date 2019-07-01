using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
#=
#1) Lar.depth_first_search:  from Tarjan, '72 (SIAM J. Computing, 1)
#2) Lar.edge_biconnect:  	from Tarjan, '72	(SIAM J. Computing, 1)
#3) Lar.biconnectedComponent:  from Hopcroft, Tarjan, '73 (ACM Communications )
=#
using ViewerGL
GL = ViewerGL

# Example 1 (DFS)
# Graph description as set of edges via EV array
EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],[6,7],[4,7],[2,7],[5,8],[1,8],[3,8]]
V = hcat([[2.,1],[2,2],[1,2],[0,3],[0,0],[3,0],[3,3],[1,1]]...)
VV = [[k] for k=1:size(V,2)]

# spanning tree: root defaults to node 1
spanningtree, fronds = Lar.depth_first_search(EV)

mesh1 = GL.numbering(1.)((V, (VV,spanningtree)), GL.COLORS[2])
mesh2 = GL.numbering(1.)((V, (VV,fronds)), GL.COLORS[4])
meshes = cat([mesh1,mesh2])
GL.VIEW( meshes );


# Example 2 (edge_biconnect vs biconnectedComponent)

V = hcat([[1.,2],[1,3],[0,3],[1,4],[2,3],[2,2],[1,1],[1,0],[0,1]]...)
EV = [[3,4],[2,4],[2,3],[2,5],[1,2],[5,6],[1,6],[1,7],[7,9],[8,9],[7,8]]
VV = [[k] for k=1:size(V,2)]

EVs = Lar.edge_biconnect(EV::Lar.Cells) # 2-edge-connected components
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );

V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );

copEV = Lar.characteristicMatrix(EV)
FE = Lar.Arrangement.biconnected_components(copEV) # 2-connected components (H & T)
EVs = [[EV[e] for e in F] for F in FE]
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );


# Example 3 (edge_biconnect vs biconnectedComponent vs biconnected_components)

V = hcat([[0,4],[4,4],[4,0],[2,0],[0,0],[2,2],[1,2],[3,2],[1,3],[3,3]]...)
EV = [[1,2],[9,10],[1,5],[7,9],[8,10],[2,3],[6,7],[6,8],[4,6],[4,5],[3,4]]
VV = [[k] for k=1:size(V,2)]

EVs = Lar.edge_biconnect(EV::Lar.Cells) # 2-edge-connected components (T)
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );

V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components (H & T)
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );

copEV = Lar.characteristicMatrix(EV)
FE = Lar.Arrangement.biconnected_components(copEV) # 2-connected components (H & T)
EVs = [[EV[e] for e in F] for F in FE]
meshes = [ GL.numbering(1.)((V, (VV,EVs[k])), GL.COLORS[k]) for k=1:length(EVs) ]
GL.VIEW( cat(meshes) );
