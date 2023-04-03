# test of unit cube with 7 faces, 15 edges, and 10 vertices.

V,(VV,EV,FV) = Lar.cuboidGrid([1,1,1],true)
W = [V [-.000000000000000001,.5,0] [-.000000000000000001,.5,1]]
W = [V [0,.5,0] [0,.5,1]]

EW = [EV; [[9,10],[2,10],[1,9]] ]
EW =[[1,2],[3,4],[5,6],[7,8],[3,9],[4,10],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8],[9,10],[2,10],[1,9]]

WW = [[k] for k=1:size(W,2)]

FW = [
    [1, 2, 9, 10],
    [3, 4, 9, 10],
    [5, 6, 7, 8],
    [1, 2, 5, 6],
    [3, 4, 7, 8],
    [1, 3, 5, 7, 9],
    [2, 4, 6, 8, 10]
]

length(WW) - length(EW) + length(FW) == 2 == (1+1) # one interior one exterior cells; 
length(WW) - length(EW) + length(FW) - 2 == 0 # Euler characteristic in E^d, d=3

GL.VIEW( GL.numbering(.5)((W,[WW, EW, FW]),GL.COLORS[1]) );
W,bases,coboundaries = Lar.chaincomplex(W,FW,EW)

julia> W
3×10 Matrix{Float64}:
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  -1.0e-18  -1.0e-18
 0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0   0.5       0.5
 0.0  1.0  0.0  1.0  0.0  1.0  0.0  1.0   0.0       1.0

julia> coboundaries[1]
15×10 SparseMatrixCSC{Int8, Int64} with 30 stored entries:
 -1   1   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
  ⋅   ⋅  -1   1   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
  ⋅  -1   ⋅   1   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
 -1   ⋅   1   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅  -1   1   ⋅   ⋅   ⋅  ⋅
  ⋅   ⋅  -1   ⋅   1   ⋅   ⋅   ⋅   ⋅  ⋅
  ⋅   ⋅   ⋅  -1   ⋅   1   ⋅   ⋅   ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  -1   1   ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  -1  1
  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  -1   ⋅   1  ⋅
  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  -1   ⋅  1
 -1   ⋅   ⋅   ⋅   ⋅   ⋅   1   ⋅   ⋅  ⋅
  ⋅  -1   ⋅   ⋅   ⋅   ⋅   ⋅   1   ⋅  ⋅
  ⋅   ⋅   ⋅   ⋅  -1   ⋅   ⋅   ⋅   1  ⋅
  ⋅   ⋅   ⋅   ⋅   ⋅  -1   ⋅   ⋅   ⋅  1

julia> coboundaries[2]
7×15 SparseMatrixCSC{Int8, Int64} with 30 stored entries:
 1  -1  1  -1  ⋅  ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
 ⋅  -1  ⋅   ⋅  1  1  -1   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅   ⋅  ⋅
 ⋅   ⋅  ⋅   ⋅  ⋅  ⋅   ⋅   1  -1  -1   1   ⋅   ⋅   ⋅  ⋅
 1   ⋅  ⋅   ⋅  ⋅  ⋅   ⋅  -1   ⋅   ⋅   ⋅  -1   1   ⋅  ⋅
 ⋅   ⋅  ⋅   ⋅  1  ⋅   ⋅   ⋅  -1   ⋅   ⋅   ⋅   ⋅  -1  1
 ⋅   ⋅  ⋅   1  ⋅  1   ⋅   ⋅   ⋅  -1   ⋅  -1   ⋅   1  ⋅
 ⋅   ⋅  1   ⋅  ⋅  ⋅   1   ⋅   ⋅   ⋅  -1   ⋅  -1   ⋅  1

julia> coboundaries[3]
2×7 SparseMatrixCSC{Int8, Int64} with 14 stored entries:
 -1   1   1   1  -1  -1   1
  1  -1  -1  -1   1   1  -1

