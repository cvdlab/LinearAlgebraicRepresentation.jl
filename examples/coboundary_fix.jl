using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

# Compute coboundary_1 in 2D via product FV * EV^t with fixing of redundancies

FV = [[1,2,3,4,5,17,16,12],
[1,2,3,4,6,7,8,9,10,11,12,13,14,15],
[4,5,9,11,12,13,14,15,16,17],
[2,3,6,7], [8,9,10,11]]

EV = [[1,2],[2,3],[3,4],[4,5],[1,12],[2,6],[3,7],[4,9],[5,17],[6,7],[8,9],
[8,10],[9,11],[10,11],[11,15],[12,13],[12,16],[13,14],[14,15],[16,17]]

V = [0   2   5   7  10   2   5   3   7  3  7  0  3  3  7  0  10;
    16  16  16  16  16  13  13  11  11  8  8  5  5  2  2  0   0]

# non-convex cells (convex=false), including outer cell (exterior=true)
copFE = Lar.coboundary_1( V::Lar.Points, FV::Lar.Cells, EV::Lar.Cells, false,true );
Matrix(copFE)

VV = [[k] for k in 1:size(V,2)];
using Plasm
Plasm.view( Plasm.numbering(3)((V,[VV, EV, FV])) )
