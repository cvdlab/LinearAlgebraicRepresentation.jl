

"""
	extrudeSimplicial(model::LAR, pattern::Array)::LAR
	
Algorithm for multimensional extrusion of a simplicial complex. 
Can be applied to 0-, 1-, 2-, ... simplicial models, to get a 1-, 2-, 3-, .... model.
The pattern `Array` is used to specify how to decompose the added dimension.

A `model` is a LAR model, i.e. a pair (vertices,cells) to be extruded, whereas pattern is an array of `Int64`, to be used as lateral measures of the *extruded* model. `pattern` elements are assumed as either *solid* or *empty* measures, according to their (+/-) sign.

# Example
```julia
julia> V = [[0,0] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];

julia> FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];

julia> pattern = repeat([1,2,-3],outer=4);

julia> model = (V,FV);

julia> W,FW = extrudeSimplicial(model, pattern);

julia> LARVIEW.view(W,FW)
```
"""
function extrudeSimplicial(model::Lar.LAR, pattern)
	V = [model[1][:,k] for k=1:size(model[1],2)]
    FV = model[2]
    d, m = length(FV[1]), length(pattern)
    coords = collect(cumsum(append!([0], abs.(pattern))))
    offset, outcells, rangelimit = length(V), [], d*m
    for cell in FV  
        tube = [v+k*offset for k in range(0, m+1) for v in cell]
        cellTube = [tube[k:k+d] for k in range(1, rangelimit)]
        outcells = vcat(outcells, reshape(cellTube, d, m))
    end
    cellGroups = []
    for i in 1:size(outcells, 2)
        if pattern[i]>0
            cellGroups = vcat(cellGroups, outcells[:, i])
        end
    end
    outVertices = [vcat(v, [z]) for z in coords for v in V]
    cellGroups = convert(Array{Array{Int, 1}, 1}, cellGroups) 
    outModel = outVertices, cellGroups
    hcat(outVertices...), cellGroups
end
function extrudeSimplicial(model::Union{Any,Lar.Cells}, pattern)
	V,FV = model
    d, m = length(FV[1]), length(pattern)
    coords = collect(cumsum(append!([0], abs.(pattern))))
    offset, outcells, rangelimit = length(V), [], d*m
    for cell in FV  
        tube = [v+k*offset for k in range(0, m+1) for v in cell]
        cellTube = [tube[k:k+d] for k in range(1, rangelimit)]
        outcells = vcat(outcells, reshape(cellTube, d, m))
    end
    cellGroups = []
    for i in 1:size(outcells, 2)
        if pattern[i]>0
            cellGroups = vcat(cellGroups, outcells[:, i])
        end
    end
    outVertices = [vcat(v, [z]) for z in coords for v in V]
    cellGroups = convert(Array{Array{Int, 1}, 1}, cellGroups) 
    outModel = outVertices, cellGroups
    hcat(outVertices...), cellGroups
end



"""
	simplexGrid(shape::Array)::LAR
	
Generate a simplicial complex decomposition of a cubical grid of ``d``-cuboids, where ``d`` is the length of `shape` array. Vertices (0-cells) of the grid have `Int64` coordinates.

# Examples
```julia
julia> simplexGrid([0]) # 0-dimensional complex
# output
([0], Array{Int64,1}[])

julia> V,EV = simplexGrid([1]) # 1-dimensional complex
# output
([0 1], Array{Int64,1}[[1, 2]])

julia> V,FV = simplexGrid([1,1]) # 2-dimensional complex
# output
([0 1 0 1; 0 0 1 1], Array{Int64,1}[[1, 2, 3], [2, 3, 4]])

julia> V,CV = simplexGrid([10,10,1]) # 3-dimensional complex
# output
([0 1 … 9 10; 0 0 … 10 10; 0 0 … 1 1], Array{Int64,1}[[1, 2, 12, 122], [2, 12, 122, 123], [12, 122, 123, 133], [2, 12, 13, 123], [12, 13, 123, 133], [13, 123, 133, 134], [2, 3, 13, 123], [3, 13, 123, 124], [13, 123, 124, 134], [3, 13, 14, 124]  …  [119, 229, 230, 240], [109, 119, 120, 230], [119, 120, 230, 240], [120, 230, 240, 241], [109, 110, 120, 230], [110, 120, 230, 231], [120, 230, 231, 241], [110, 120, 121, 231], [120, 121, 231, 241], [121, 231, 241, 242]])

julia> V
# output
3×242 Array{Int64,2}:
 0  1  2  3  4  5  6  7  8  9  10  0  1  2  3  …   1   2   3   4   5   6   7   8   9  10
 0  0  0  0  0  0  0  0  0  0   0  1  1  1  1     10  10  10  10  10  10  10  10  10  10
 0  0  0  0  0  0  0  0  0  0   0  0  0  0  0      1   1   1   1   1   1   1   1   1   1

julia> using LARVIEW

julia> hpc = LARVIEW.lar2exploded_hpc(V,CV) # exploded visualization of the grid

julia> LARVIEW.view(hpc)

julia> V,HV = simplexGrid([1,1,1,1]) # 4-dim cellular complex from the 4D simplex
# output
([0 1 … 0 1; 0 0 … 1 1; 0 0 … 1 1; 0 0 … 1 1], Array{Int64,1}[[1, 2, 3, 5, 9], [2, 3, 5, 9, 10], [3, 5, 9, 10, 11], [5, 9, 10, 11, 13], [2, 3, 5, 6, 10], [3, 5, 6, 10, 11], [5, 6, 10, 11, 13], [6, 10, 11, 13, 14], [3, 5, 6, 7, 11], [5, 6, 7, 11, 13]  …  [4, 6, 10, 11, 12], [6, 10, 11, 12, 14], [3, 4, 6, 7, 11], [4, 6, 7, 11, 12], [6, 7, 11, 12, 14], [7, 11, 12, 14, 15], [4, 6, 7, 8, 12], [6, 7, 8, 12, 14], [7, 8, 12, 14, 15], [8, 12, 14, 15, 16]])
```
"""
function simplexGrid(shape)
    model = [[]], [[1]]
    for item in shape
        model = extrudeSimplicial(model, fill(1, item))
    end
    V, CV = model
    V = convert(Array{Float64,2}, V)
    return V, CV
end





"""
	simplexFacets(simplices::Cells)::Cells

Compute the `(d-1)`-skeleton (set of `facets`) of a simplicial `d`-complex.

# Examples
```julia
julia> V,FV = simplexGrid([1,1]) # 2-dimensional complex
# output
([0 1 0 1; 0 0 1 1], Array{Int64,1}[[1, 2, 3], [2, 3, 4]])

julia> LARVIEW.view(V,FV)

julia> W,CW = extrudeSimplicial((V,FV), [1])
([0.0 1.0 … 0.0 1.0; 0.0 0.0 … 1.0 1.0; 0.0 0.0 … 1.0 1.0], 
Array{Int64,1}[[1,2,3,5],[2,3,5,6],[3,5,6,7],[2,3,4,6],[3,4,6,7],[4,6,7,8]])

julia> FW = simplexFacets(CW)
18-element Array{Any,1}:
[[1,3,5],[5,6,7],[3,5,7],[3,6,7],[4,6,7],[4,7,8],[4,6,8],
[6,7,8],[3,5,6],[2,3,5],[2,3,4],[3,4,7],[1,2,3],[2,4,6],[2,5,6],
[1,2,5],[2,3,6],[3,4,6]]

julia> LARVIEW.view(W,FW)
```
"""
function simplexFacets(simplices)
    out = []
    d = length(simplices[1])
    for simplex in simplices
    	non_oriented_facets = [ append!(simplex[1:k-1],simplex[k+1:d]) for k=1:d ]
    	if (-1)^d == 1
    		oriented_facets = non_oriented_facets
    	else
    		oriented_facets = [ length(f)==2 ? [f[2],f[1]] : [f[2],f[1],f[3:end]]  
    			for f in non_oriented_facets ]
    	end
        append!(out, oriented_facets)
    end
    out = map(sort, collect(Set(out)))
    out = sort(out, by=x->(x[1],x[2]))
    return convert(Lar.Cells, out)
end


