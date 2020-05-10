### Two cubes arrangement

Here we first create the `Model` of a cube, actually a wrapping of vertices, edges and faces of a cellular complex created by the `Lar.cuboid` primitive, with minimum point `[-1,-1,-1]`, maximum point [1,1,1], and `true` flag for generation of all bases of its boundary complex.

```
julia> V,(VV,EV,FV,CV) = Lar.cuboid([1,1,1],true,[-1,-1,-1])
```

Then we generate the `cube::Model` using the easiest generating method :

```
julia> cube = CAGD.Model(V,[EV,FV]);

```
The `cube` *geometry* in $E^3$ is given by the $3\times 8$ `Float64` matrix
```
julia> cube.G
3×8 Array{Float64,2}:
 -1.0  -1.0  -1.0  -1.0   1.0   1.0   1.0  1.0
 -1.0  -1.0   1.0   1.0  -1.0  -1.0   1.0  1.0
 -1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0  1.0
```
The internal representation of the `cube.T` *topology* is stored as:

```
julia> typeof(cube.T)
Array{SparseMatrixCSC{Int8,Int64},1}

julia> @show map(SparseArrays.findnz, cube.T);
map(SparseArrays.findnz, cube.T) = Tuple{Array{Int64,1},Array{Int64,1},Array{Int8,1}}[([1, 5, 9, 1, 6, 10, 2, 5, 11, 2, 6, 12, 3, 7, 9, 3, 8, 10, 4, 7, 11, 4, 8, 12], [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8], [-1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1]), ([1, 3, 1, 4, 2, 3, 2, 4, 1, 5, 1, 6, 2, 5, 2, 6, 3, 5, 3, 6, 4, 5, 4, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12], [1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1]), ([], [], [])]
```
We can see the `Matrix` representation:

```
julia> map(Matrix, cube.T)
3-element Array{Array{Int8,2},1}:
 [-1 1 … 0 0; 0 0 … 0 0; … ; 0 0 … 1 0; 0 0 … 0 1]
 [1 -1 … 0 0; 0 0 … 0 0; … ; 0 0 … 1 0; 0 0 … 0 1]
 Array{Int8}(undef,0,6)                           
```
Or, more explicitly: $[\partial_0]$
```
julia> Matrix(cube.T[1])
12×8 Array{Int8,2}:
 -1   1   0   0   0   0   0  0
  0   0  -1   1   0   0   0  0
  0   0   0   0  -1   1   0  0
  0   0   0   0   0   0  -1  1
 -1   0   1   0   0   0   0  0
  0  -1   0   1   0   0   0  0
  0   0   0   0  -1   0   1  0
  0   0   0   0   0  -1   0  1
 -1   0   0   0   1   0   0  0
  0  -1   0   0   0   1   0  0
  0   0  -1   0   0   0   1  0
  0   0   0  -1   0   0   0  1
```
and $[\partial_1]$
```
julia> Matrix(cube.T[2])
6×12 Array{Int8,2}:
 1  -1   0   0  -1  1   0   0   0   0   0  0
 0   0   1  -1   0  0  -1   1   0   0   0  0
 1   0  -1   0   0  0   0   0  -1   1   0  0
 0   1   0  -1   0  0   0   0   0   0  -1  1
 0   0   0   0   1  0  -1   0  -1   0   1  0
 0   0   0   0   0  1   0  -1   0  -1   0  1
```
