using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

# Compute coboundary_1 in 2D via product FV * EV^t with fixing of redundancies.
# The cellular complex includes non-convex cells.
# The coboundary_1 matrix is generated coherently oriented, and with the outer cell on the first row.

V = [5. 0  6  9  5  3  0  5  8  10  10  10  0  0  3  9  8  6  5;
     9  9  5  5  5  6  6  0  0   9   5   0  0  3  3  8  3  8  3]

EV = [[1,2],[3,4],[3,5],[1,10],[16,18],[7,14],[6,15],[8,13],[13,14],[10,11],
[11,12],[8,19],[9,17],[1,5],[17,19],[6,7],[14,15],[4,11],[3,18],[4,16],[8,9],[9,12],[2,7]];

FV = [[1,2,7,6,15,14,13,8,19,17,9,12,11,4,3,5],
[3,4,16,18],[6,7,14,15],[1,3,4,5,10,11,16,18],
[8,9,17,19],[1,2,7,14,13,8,9,12,11,10]]

VV = [[k] for k in 1:size(V,2)];
GL.VIEW( GL.numbering(3)((V,[VV, EV, FV])) )

# testing Lar.coboundary_1  #TODO DEBUG
convex = false;
exterior = true;
copFE = Lar.coboundary_1( V::Lar.Points, FV::Lar.Cells, EV::Lar.Cells, convex,exterior );

[ collect(1:size(copFE,2))'; Matrix(copFE) ]
