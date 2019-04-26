using LinearAlgebraicRepresentation
using Plasm
using SparseArrays
Lar = LinearAlgebraicRepresentation
Lara = Lar.Arrangement


V = [
	0.0 0.4 0.5 0.6 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.4 0.5 0.6 1.0;
	0.0 0.0 0.0 0.0 0.0 0.4 0.4 0.5 0.5 0.6 0.6 1.0 1.0 1.0 1.0 1.0
]

EV = [
	[ 1,  2], [ 1,  6], [ 2,  6],
	[ 4,  5], [ 4,  7], [ 5,  7],
	[10, 12], [10, 13], [12, 13],
	[11, 15], [11, 16], [15, 16],
	[ 3,  8], [ 3,  9], [ 3, 14], [ 8,  9], [ 8, 14], [ 9, 14]
]

## 1 - Input Reading
Plasm.view(Plasm.numbering(0.5)((V,[[[k] for k=1:size(V,2)], EV])));


# Planar Arrangement 1
W = convert(Lar.Points, V'); # Infering type for W = V'
copEV = Lar.coboundary_0(EV::Lar.Cells);
W1, copEV1 = Lar.planar_arrangement_1(W::Lar.Points, copEV::Lar.ChainOp);


## 2 - 2-cells fragmentation
EV1 = Lar.cop2lar(copEV1);
V1 = convert(Lar.Points, W1');
Plasm.view(Plasm.numbering(0.5)((V1,[[[k] for k=1:size(V1,2)], EV1])));


# Biconnected COmponent Evaluation
bicon_comps = Lar.Arrangement.biconnected_components(copEV1);


## 3 - Biconnected Components
hpcs = [ Plasm.lar2hpc(V1,[EV1[e] for e in comp]) for comp in bicon_comps ]
Plasm.view([
	Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k])
	for k = 1 : (length(hpcs))
])


# computation of 2-cells and 2-boundary

n = size(bicon_comps, 1)
shells = Array{Lar.Chain, 1}(undef, n)
boundaries = Array{Lar.ChainOp, 1}(undef, n)
EVs = Array{Lar.ChainOp, 1}(undef, n)
# for each component
for p=1:n
	ev = copEV1[sort(bicon_comps[p]), :]
	fe = Lara.minimal_2cycles(W1, ev)
	global shell_num = Lara.get_external_cycle(W1, ev, fe)
	EVs[p] = ev
	global tokeep = setdiff(1:fe.m, shell_num)
	boundaries[p] = fe[tokeep, :]
	shells[p] = fe[shell_num, :]
end

shell_bboxes = []
for i in 1:n
	vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
	push!(shell_bboxes, Lar.bbox(W1[vs_indexes, :]))
end

containment_graph = Lara.pre_containment_test(shell_bboxes)
containment_graph = Lara.prune_containment_graph(
	n, W1, EVs, shells, containment_graph
)
Lara.transitive_reduction!(containment_graph)
W2 = W1;
copEV2, copFE2 = Lara.cell_merging(
	n, containment_graph, W2, EVs, boundaries, shells, shell_bboxes
)


## 4 - 3-cells identification & dangling 1-cells elimination
Plasm.view( Plasm.numbering1(0.5)((W2, copEV2, copFE2)) )


# 5 - Colorfull Representation
triangulated_faces = Lar.triangulate2D(W2, [copEV2, copFE2])
V2 = convert(Lar.Points, W2')
FVs2 = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewcolor(V2::Lar.Points, FVs2::Array{Lar.Cells})

# polygonal face boundaries
EVs2 = Lar.FV2EVs(copEV2, copFE2)
EVs2 = convert(Array{Array{Array{Int64,1},1},1}, EVs2)
Plasm.viewcolor(V2::Lar.Points, EVs2::Array{Lar.Cells})

# 6 - Exploded Representation
model = V2,EVs2
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
