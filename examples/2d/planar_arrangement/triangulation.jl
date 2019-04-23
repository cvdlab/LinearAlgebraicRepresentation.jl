using LinearAlgebraicRepresentation
using Plasm
using SparseArrays
Lar = LinearAlgebraicRepresentation
Lara = Lar.Arrangement

# Data Reading
V = [
	0.0 0.5 3.0 5.0 2.0 2.0 0.0 1.0 0.0 3.0; 
	1.0 1.0 1.0 1.0 2.0 0.0 3.0 1.5 0.0 1.5
	];
EV = [[1, 4], [2, 3], [5, 6], [2, 5], [3, 6], [7, 10], [8, 10], [9, 10]];

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
W2, copEV2, copFE2 = Lar.planar_arrangement_2(W1, copEV1, bicon_comps)

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