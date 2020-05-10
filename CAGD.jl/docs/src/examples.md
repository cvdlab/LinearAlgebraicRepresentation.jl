## Cube generation as a `Model`

Here we first create the `Model` of a cube, actually a wrapping of vertices, edges and faces of a cellular complex created by the `Lar.cuboid` primitive, with minimum point `[-1,-1,-1]`, maximum point [1,1,1], and `true` flag for generation of all bases of its boundary complex.

```julia
julia> V,(VV,EV,FV,CV) = Lar.cuboid([1,1,1],true,[-1,-1,-1])
```

Then we generate the `cube::Model` using the easiest generating method :

```julia
julia> cube = CAGD.Model(V,[EV,FV]);
```

The `cube` *geometry* in ``E^3`` is given by the $3\times 8$ `Float64` matrix

```julia
julia> cube.G
38 Array{Float64,2}:
 -1.0  -1.0  -1.0  -1.0   1.0   1.0   1.0  1.0
 -1.0  -1.0   1.0   1.0  -1.0  -1.0   1.0  1.0
 -1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0  1.0
```

The internal representation of the `cube.T` *topology* is stored as:

```julia
julia> typeof(cube.T)
Array{SparseMatrixCSC{Int8,Int64},1}

julia> @show map(SparseArrays.findnz, cube.T);
map(SparseArrays.findnz, cube.T) = Tuple{Array{Int64,1},Array{Int64,1},Array{Int8,1}}[([1, 5, 9, 1, 6, 10, 2, 5, 11, 2, 6, 12, 3, 7, 9, 3, 8, 10, 4, 7, 11, 4, 8, 12], [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 8, 8], [-1, -1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1]), ([1, 3, 1, 4, 2, 3, 2, 4, 1, 5, 1, 6, 2, 5, 2, 6, 3, 5, 3, 6, 4, 5, 4, 6], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12], [1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1]), ([], [], [])]
```
We can see the `Matrix` representation:

```julia
julia> map(Matrix, cube.T)
3-element Array{Array{Int8,2},1}:
 [-1 1  0 0; 0 0  0 0;  ; 0 0  1 0; 0 0  0 1]
 [1 -1  0 0; 0 0  0 0;  ; 0 0  1 0; 0 0  0 1]
 Array{Int8}(undef,0,6)                           
```

Or, more explicitly: $[\delta_0] : `V` \to `E`$

```julia
julia> Matrix(cube.T[1])
128 Array{Int8,2}:
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

and $[\delta_1] : `E` \to `F`$

```julia
julia> Matrix(cube.T[2])
612 Array{Int8,2}:
 1  -1   0   0  -1  1   0   0   0   0   0  0
 0   0   1  -1   0  0  -1   1   0   0   0  0
 1   0  -1   0   0  0   0   0  -1   1   0  0
 0   1   0  -1   0  0   0   0   0   0  -1  1
 0   0   0   0   1  0  -1   0  -1   0   1  0
 0   0   0   0   0  1   0  -1   0  -1   0  1
```

## `Struct` and `Model` assembly of two tuple instances.  

- `Struct` type is used to aggregate hierarchically in world coordinates both geometric tuples and other `Struct` objects, defined in local coordinates, and transformed eventually by affine maps, using homogeneous coordinates.
- `Model` type is used to wrap geometry and topology of an object together, and jointly assign to a variable symbol.

In a next release `Struct` will be able to combine hierarchically any number of `Model` objects.

```julia
julia> assembly = Lar.Struct([ (V,EV,FV), Lar.t(1,1,1), (V,EV,FV) ]);

julia> typeof(assembly)
LinearAlgebraicRepresentation.Struct

julia> V,EV,FV = Lar.struct2lar(assembly);

julia> twocubes = CAGD.Model(V, [EV,FV]);
```

where we have, of course, `8+8` vertices, `12+12` edges, and `6+6` faces.  Notice that the number of 3-cell is yet undetermined (0), because `twocubes` is not a cellular complex (missing the boundary compatibility) but a collection of two (intersecting) cellular complexes. 

```julia
julia> twocubes.G
316 Array{Float64,2}:
 -1.0  -1.0  -1.0  -1.0   1.0   1.0    0.0  0.0  0.0  2.0  2.0  2.0  2.0
 -1.0  -1.0   1.0   1.0  -1.0  -1.0     0.0  2.0  2.0  0.0  0.0  2.0  2.0
 -1.0   1.0  -1.0   1.0  -1.0   1.0     2.0  0.0  2.0  0.0  2.0  0.0  2.0
 
julia> map(size, twocubes.T)
3-element Array{Tuple{Int64,Int64},1}:
 (24, 16)
 (12, 24)
 (0, 12) 
```
In order to create a well defined cellular complex, we must apply a pipeline of transformations, starting from object's fragmentation. Every face is split in 2D against all faces of possible intersection with it (using interval-trees). Solid cells are not yet reconstructed:

```julia
julia> split_model = CAGD.facesplitting(twocubes)
CAGD.Model([-1.0 -1.0 … 2.0 2.0; -1.0 -1.0 … 0.0 2.0; -1.0 1.0 … 2.0 2.0], SparseArrays.SparseMatrixCSC{Int8,Int64}[ ....

julia> size(split_model.G)
(3, 66)

julia> map(size,split_model.T)
3-element Array{Tuple{Int64,Int64},1}:
 (72, 66)
 (18, 72)
 (0, 18) 
```
Note that numbers of vertices, edges, faces, and cells were now 66, 72, 18, 0, respectively, since many are computed *independently* in multiple instances, starting from each face in 2D:

```julia
julia> congruent_model = CAGD.mergemodel(split_model, signed_merge=true)

julia> displayModel(congruent_model)

julia> size(congruent_model.G)
(3, 22)

julia> map(size, congruent_model.T)
3-element Array{Tuple{Int64,Int64},1}:
 (36, 22)
 (18, 36)
 (0, 18) 
```
The numbers of vertices, edges, faces, and cells were now 2, 36, 18, 0, after having identified ``\epsilon``-congruent vertices, edges, and faces.
