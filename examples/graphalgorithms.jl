using LinearAlgebraicRepresentation
using Plasm, DataStructures
using PyCall
Lar = LinearAlgebraicRepresentation
p = PyCall.pyimport("pyplasm")

# Example 1

EV = [[1,2],[2,3],[3,4],[4,5],[5,6],[1,6],
[6,7],[4,7],[2,7],[5,8],[1,8],[3,8]]

VV = Lar.verts2verts(EV)
V = hcat([[2.,1],[2,2],[1,2],[0,3],[0,0],[3,0],[3,3],[1,1]]...)

spanningtree, fronds = Lar.depth_first_search(EV)

hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])

# Example 2

spanningtree, fronds = Lar.depth_first_search(EV,5)

hpc1 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], spanningtree]))
hpc2 = Plasm.numbering(1.25)((V,[[[k] for k=1:size(V,2)], fronds]))
Plasm.view([ Plasm.color("cyan")(hpc1) , Plasm.color("magenta")(hpc2) ])


