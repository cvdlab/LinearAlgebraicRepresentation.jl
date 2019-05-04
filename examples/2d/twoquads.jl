using LinearAlgebraicRepresentation
using Plasm, SparseArrays
Lar = LinearAlgebraicRepresentation

v = [-1.0 0.0 0.0 1.0 -1.0 0.0 0.0 1.0 -1.0 0.0 -1.0 0.0 1.0 0.0 0.0 1.0; 1.0 0.0 1.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.0 0.0 0.0 1.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16]]
Plasm.view(Plasm.numbering(0.5)((v,[[[k] for k=1:size(v,2)],ev])))

verts = convert(Lar.Points,v')
copEV = Lar.coboundary_0(convert(Lar.Cells,ev))

V, copEV, copFE = Lar.Arrangement.planar_arrangement(verts::Lar.Points, copEV::Lar.ChainOp)

W = convert(Lar.Points,V')
EW = Lar.cop2lar(copEV)
FE = Lar.cop2lar(copFE)
FW = [collect(Set(vcat(cat([EW[f] for e in f])...))) for f in FE]
Plasm.view(Plasm.numbering(0.5)((W,[[[k] for k=1:size(W,2)],EW, FW])))
