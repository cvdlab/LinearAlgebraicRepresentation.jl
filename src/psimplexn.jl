#####################################################
### From Julia 0.6 to Julia 1.0
### --- SIMPLEXN - parallel version ---
### Authors: Fabio Fatelli, Emanuele Loprevite
#####################################################

"""
	plarExtrude1(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}, pattern::Array{Int64,1}) where T<:Real
		::Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}

Algorithm for multidimensional extrusion of a simplicial complex. 
Can be applied to 0-, 1-, 2-, ... simplicial models, to get a 1-, 2-, 3-, ... model.
The pattern `Array` is used to specify how to decompose the added dimension.

A `model` is a LAR model, i.e. a pair (vertices,cells) to be extruded, whereas pattern is an array of `Int64`, to be used as lateral measures of the *extruded* model. `pattern` elements are assumed as either *solid* or *empty* measures, according to their (+/-) sign.

# Example
```julia
julia> V = [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1], [0,2], [1,2], [2,2]];

julia> FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];

julia> pattern = repeat([1,2,-3],outer=4);

julia> model = (V,FV);

julia> W,FW = plarExtrude1(model, pattern);
```
"""
# Generation of the output model vertices in a multiple extrusion of a LAR model
using SharedArrays

@everywhere function plarExtrude1(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}, pattern::Array{Int64,1}) where T<:Real
	V, FV = model
	d, m = length(FV[1]), length(pattern)
	dd1 = d*(d+1)
	coords = cumsum(append!([0],abs.(pattern))) # built-in function cumsum
	outVertices = @spawn [vcat(v,z) for z in coords for v in V]
	offset, outcells, rangelimit = length(V), SharedArray{Int64}(m,dd1*length(FV)), d*m
	@sync @distributed for j in 1:length(FV)
		tube = [v+k*offset for k in 0:m for v in FV[j]]
		celltube = Int64[]
		celltube = @sync @distributed (append!) for k in 1:rangelimit
			tube[k:k+d]
		end
		outcells[:,(j-1)*dd1+1:(j-1)*dd1+dd1] = reshape(celltube,dd1,m)'
	end
	p = convert(SharedArray,findall(x->x>0,pattern))
	cellGroups = SharedArray{Int64}(length(p),size(outcells)[2])
	@sync @distributed for k in 1:length(p)
		cellGroups[k,:] = outcells[p[k],:]
	end
	cellGroupsL = vec(cellGroups')
	outCellGroups = Array{Int64,1}[]
	outCellGroups = @distributed (append!) for k in 1:d+1:length(cellGroupsL)
			[cellGroupsL[k:k+d]]
		end
	return fetch(outVertices), outCellGroups
end

"""
	plarSimplexGrid1(shape::Array{Int64,1})::Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}

Generate a simplicial complex decomposition of a cubical grid of ``d``-cuboids, where ``d`` is the length of `shape` array. Vertices (0-cells) of the grid have `Int64` coordinates.

# Example
```julia
julia> plarSimplexGrid1([0]) # 0-dimensional complex
# output
(Array{Int64,1}[[0]], Array{Int64,1}[])

julia> V,EV = plarSimplexGrid1([1]) # 1-dimensional complex
# output
(Array{Int64,1}[[0], [1]], Array{Int64,1}[[0, 1]])

julia> V,FV = plarSimplexGrid1([1,1]) # 2-dimensional complex
# output
(Array{Int64,1}[[0, 0], [1, 0], [0, 1], [1, 1]], Array{Int64,1}[[0, 1, 2], [1, 2, 3]])

julia> V,CV = plarSimplexGrid1([10,10,1]) # 3-dimensional complex
# output
(Array{Int64,1}[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0], [5, 0, 0], [6, 0, 0], [7, 0, 0], [8, 0, 0], [9, 0, 0]  …  [1, 10, 1], [2, 10, 1], [3, 10, 1], [4, 10, 1], [5, 10, 1], [6, 10, 1], [7, 10, 1], [8, 10, 1], [9, 10, 1], [10, 10, 1]], Array{Int64,1}[[0, 1, 11, 121], [1, 11, 121, 122], [11, 121, 122, 132], [1, 11, 12, 122], [11, 12, 122, 132], [12, 122, 132, 133], [1, 2, 12, 122], [2, 12, 122, 123], [12, 122, 123, 133], [2, 12, 13, 123]  …  [118, 228, 229, 239], [108, 118, 119, 229], [118, 119, 229, 239], [119, 229, 239, 240], [108, 109, 119, 229], [109, 119, 229, 230], [119, 229, 230, 240], [109, 119, 120, 230], [119, 120, 230, 240], [120, 230, 240, 241]])

julia> V,HV = plarSimplexGrid1([1,1,1,1]) # 4-dim cellular complex from the 4D simplex
# output
(Array{Int64,1}[[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0], [0, 1, 1, 0], [1, 1, 1, 0], [0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 0, 1], [0, 0, 1, 1], [1, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1]], Array{Int64,1}[[0, 1, 2, 4, 8], [1, 2, 4, 8, 9], [2, 4, 8, 9, 10], [4, 8, 9, 10, 12], [1, 2, 4, 5, 9], [2, 4, 5, 9, 10], [4, 5, 9, 10, 12], [5, 9, 10, 12, 13], [2, 4, 5, 6, 10], [4, 5, 6, 10, 12]  …  [3, 5, 9, 10, 11], [5, 9, 10, 11, 13], [2, 3, 5, 6, 10], [3, 5, 6, 10, 11], [5, 6, 10, 11, 13], [6, 10, 11, 13, 14], [3, 5, 6, 7, 11], [5, 6, 7, 11, 13], [6, 7, 11, 13, 14], [7, 11, 13, 14, 15]])
```
"""
# Generation of simplicial grids of any dimension and shape
@everywhere function plarSimplexGrid1(shape::Array{Int64,1})
	model = [Int64[]],[[0]] # the empty simplicial model
	for item in shape # no parallel
		model = plarExtrude1(model,repeat([1],item))
	end
	return model
end


"""
	plarSimplexFacets(simplices::Array{Array{Int64,1},1})::Array{Array{Int64,1},1}

Compute the `(d-1)`-skeleton (set of `facets`) of a simplicial `d`-complex.

# Example
```julia
julia> V,FV = plarSimplexGrid1([1,1])
# output
(Array{Int64,1}[[0, 0], [1, 0], [0, 1], [1, 1]], Array{Int64,1}[[0, 1, 2], [1, 2, 3]])

julia> W,CW = plarExtrude1((V,FV), [1])
# output
(Array{Int64,1}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], Array{Int64,1}[[0, 1, 2, 4], [1, 2, 4, 5], [2, 4, 5, 6], [1, 2, 3, 5], [2, 3, 5, 6], [3, 5, 6, 7]])

julia> FW = plarSimplexFacets(CW)
18-element Array{Array{Int64,1},1}:
 [0, 1, 2]
 [0, 1, 4]
 [0, 2, 4]
 [1, 2, 3]
 [1, 2, 4]
 [1, 2, 5]
 [1, 3, 5]
 [1, 4, 5]
 [2, 3, 5]
 [2, 3, 6]
 [2, 4, 5]
 [2, 4, 6]
 [2, 5, 6]
 [3, 5, 6]
 [3, 5, 7]
 [3, 6, 7]
 [4, 5, 6]
 [5, 6, 7]
```
"""
# Extraction of non-oriented (d-1)-facets of d-dimensional simplices
@everywhere using Combinatorics # for combinations() function

@everywhere function plarSimplexFacets(simplices::Array{Array{Int64,1},1})
	out = Array{Int64,1}[]
	d = length(simplices[1])
	out = @distributed (append!) for simplex in simplices
			collect(combinations(simplex,d-1))
		end
	return sort!(unique(out),lt=isless)
end


"""
	pquads2tria(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}) where T<:Real
		::Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}

Give the conversion of a LAR boundary representation (B-Rep), i.e. a LAR model `V, FV` made of 2D faces, usually quads but also general polygons, into a LAR model `verts, triangles` made by triangles.

# Example
```julia
julia> model = plarSimplexGrid1([2,2,2])
# output
(Array{Int64,1}[[0, 0, 0], [1, 0, 0], [2, 0, 0], [0, 1, 0], [1, 1, 0], [2, 1, 0], [0, 2, 0], [1, 2, 0], [2, 2, 0], [0, 0, 1]  …  [2, 2, 1], [0, 0, 2], [1, 0, 2], [2, 0, 2], [0, 1, 2], [1, 1, 2], [2, 1, 2], [0, 2, 2], [1, 2, 2], [2, 2, 2]], Array{Int64,1}[[0, 1, 3, 9], [1, 3, 9, 10], [3, 9, 10, 12], [1, 3, 4, 10], [3, 4, 10, 12], [4, 10, 12, 13], [1, 2, 4, 10], [2, 4, 10, 11], [4, 10, 11, 13], [2, 4, 5, 11]  …  [15, 21, 22, 24], [13, 15, 16, 22], [15, 16, 22, 24], [16, 22, 24, 25], [13, 14, 16, 22], [14, 16, 22, 23], [16, 22, 23, 25], [14, 16, 17, 23], [16, 17, 23, 25], [17, 23, 25, 26]])

julia> pquads2tria(model)
# output
(Array{Float64,1}[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [0.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0], [0.0, 0.0, 1.0]  …  [0.25, 1.5, 1.75], [0.75, 1.5, 1.25], [0.5, 1.75, 1.5], [0.75, 1.75, 1.75], [1.25, 1.25, 1.25], [1.5, 1.25, 1.5], [1.25, 1.5, 1.75], [1.75, 1.5, 1.25], [1.5, 1.75, 1.5], [1.75, 1.75, 1.75]], Array{Int64,1}[[27, 3, 9], [27, 9, 0], [27, 0, 1], [27, 1, 3], [28, 9, 10], [28, 10, 1], [28, 1, 3], [28, 3, 9], [29, 10, 12], [29, 12, 3]  …  [72, 16, 17], [72, 17, 23], [73, 23, 25], [73, 25, 16], [73, 16, 17], [73, 17, 23], [74, 26, 25], [74, 25, 17], [74, 17, 23], [74, 23, 26]])
```
"""
# Transformation to triangles by sorting circularly the vertices of faces
using LinearAlgebra

@everywhere function pquads2tria(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}) where T<:Real
	V, FV = model
	if typeof(V) != Array{Array{Float64,1},1}
		V = convert(Array{Array{Float64,1},1},V)
	end
	out = Array{Int64,1}[]
	nverts = length(V)-1
	for face in FV # no parallel
		arr = [V[v+1] for v in face]
		centroid = sum(arr)/length(arr)
		append!(V,[centroid])
		nverts += 1
		v1, v2 = V[face[1]+1]-centroid, V[face[2]+1]-centroid
		v3 = cross(v1,v2)
		if norm(v3) < 1/(10^3)
			v1, v2 = V[face[1]+1]-centroid, V[face[3]+1]-centroid
			v3 = cross(v1,v2)
		end
		transf = inv(copy(hcat(v1,v2,v3)'))
		verts = [(V[v+1]'*transf)'[1:end-1] for v in face]
		tcentroid = sum(verts)/length(verts)
		tverts = pmap(x->x-tcentroid,verts)
		iterator = collect(zip(tverts,face))
		rverts = [[atan(reverse(iterator[i][1])...),iterator[i][2]] for i in 1:length(iterator)]
		rvertsS = sort(rverts,lt=(x,y)->isless(x[1],y[1]))
		ord = [pair[2] for pair in rvertsS]
		append!(ord,ord[1])
		edges = [[i[2],ord[i[1]+1]] for i in enumerate(ord[1:end-1])]
		triangles = pmap(x->prepend!(x,nverts),edges)
		append!(out,triangles)
	end
	return V, out
end
