using LinearAlgebraicRepresentation
using Plasm, DataStructures
using PyCall
Lar = LinearAlgebraicRepresentation
const colorkey = ["orange", "red", "green", "blue", "cyan", "magenta", "yellow", 
	"white", "black", "purple", "gray", "brown"]
#=
#1) Lar.depth_first_search:  from Tarjan, '72 (SIAM J. Computing, 1)
#2) Lar.edge_biconnect:  	from Tarjan, '72	(SIAM J. Computing, 1)
#3) Lar.biconnectedComponent:  from Hopcroft, Tarjan, '73 (ACM Communications )
=#

# Example 1 (DFS)

# Graph description as set of edges via EV array
EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],[6,7],[4,7],[2,7],[5,8],[1,8],[3,8]]
V = hcat([[2.,1],[2,2],[1,2],[0,3],[0,0],[3,0],[3,3],[1,1]]...)

spanningtree, fronds = Lar.depth_first_search(EV) # spanning tree: default to node 1
hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])

spanningtree, fronds = Lar.depth_first_search(EV,5) # spanning tree:  node 5
hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])


# Example 2 (edge_biconnect vs biconnectedComponent)

V = hcat([[1.,2],[1,3],[0,3],[1,4],[2,3],[2,2],[1,1],[1,0],[0,1]]...)
EV = [[3,4],[2,4],[2,3],[2,5],[1,2],[5,6],[1,6],[1,7],[7,9],[8,9],[7,8]]

EVs = Lar.edge_biconnect(EV::Lar.Cells) # 2-edge-connected components
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])

V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])

copEV = Lar.characteristicMatrix(EV)
FE = Lar.Arrangement.biconnected_components(copEV) # 2-connected components (H & T)
EVs = [[EV[e] for e in F] for F in FE]
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])


# Example 3 (edge_biconnect vs biconnectedComponent vs biconnected_components)

V = hcat([[0,4],[4,4],[4,0],[2,0],[0,0],[2,2],[1,2],[3,2],[1,3],[3,3]]...)
EV = [[1,2],[9,10],[1,5],[7,9],[8,10],[2,3],[6,7],[6,8],[4,6],[4,5],[3,4]]

EVs = Lar.edge_biconnect(EV::Lar.Cells) # 2-edge-connected components (T)
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])

V,EVs = Lar.biconnectedComponent((V,EV::Lar.Cells)) # 2-connected components (H & T)
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])

copEV = Lar.characteristicMatrix(EV)
FE = Lar.Arrangement.biconnected_components(copEV) # 2-connected components (H & T)
EVs = [[EV[e] for e in F] for F in FE]
hpcs = [ Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)],EVs[i]])) for i=1:length(EVs)]
Plasm.view([ Plasm.color(colorkey[k])(hpcs[k]) for k=1:length(hpcs) ])


# Example 4

V = hcat([[0.,0],[1,0],[1,1],[0,1],[2,1]]...);
EV = [[1,2],[2,3],[3,4],[4,1],[1,5]];
model = V,EV
W,EW = Lar.fragmentlines(model)
V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)
Plasm.view(Plasm.numbering()((V,[[[k] for k=1:size(V,2)],EVs[1]])))

