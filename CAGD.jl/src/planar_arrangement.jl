Lar = LinearAlgebraicRepresentation
#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1 METHODS
#    frag_edge
#    intersect_edges
#    merge_vertices!
#-------------------------------------------------------------------------------

"""
	frag_edge(
		model::CAGD.Model,
		edge_idx::Int,
		bigPI::Array{Array{Int64,1},1}
	)::Tuple{Lar.Points, Lar.ChainOp}

Splits the `edge_idx`-th edge of the Edge Topology `model.T[1]`.

This method splits the `edge_idx`-th edge between the Edges in the Topology
`EV = model.T[1]` into several parts by confronting it with the others that
intersect its bounding box (see also [`CAGD.intersect_edges`](@ref)).
Returns a set made of the detected vertices over the edge (with redundancies)
and the associated cochain (without redundancies).

See also:
 - [`CAGD.planar_arrangement_1`](@ref)
 - [`CAGD.frag_edge_channel`](@ref)
It uses: [`CAGD.intersect_edges`](@ref)
---
# Examples
```jldoctest
julia> model = CAGD.Model(hcat([
	[1.0, 0.0], [0.0, 1.0], [0.0, 0.5], [0.5, 1.0], [1.0, 1.0]
]...));
julia> CAGD.addModelCells!(model, 1, [[1, 2], [2, 5], [3, 4], [4, 5]])
julia> bigPI = CAGD.spaceIndex(model);
julia> CAGD.frag_edge(model, 1, bigPI)[1]
4×2 Array{Float64,2}:
 1.0   0.0
 0.0   1.0
 0.25  0.75
 0.0   1.0
julia> Matrix(CAGD.frag_edge(model, 1, bigPI)[2])
2×4 Array{Int8,2}:
 1  0  1  0
 0  1  1  0
```
"""
function frag_edge(
		model::CAGD.Model,
		edge_idx::Int,
		bigPI::Array{Array{Int64,1},1}
	)::Tuple{Lar.Points, Lar.ChainOp}
	V = convert(Lar.Points, model.G')
	EV = model.T[1]
	alphas = Dict{Float64, Int}()
	edge = EV[edge_idx, :]
	verts = V[edge.nzind, :]
	for i in bigPI[edge_idx]
		if i != edge_idx # && edge_idx in bigPI[i]
			intersection = CAGD.intersect_edges(V, edge, EV[i, :])
			for (point, alpha) in intersection
				verts = [verts; point]
				alphas[alpha] = size(verts, 1)
			end
		end
	end
	alphas[0.0], alphas[1.0] = [1, 2]
	alphas_keys = sort(collect(keys(alphas)))
	edge_num = length(alphas_keys)-1
	verts_num = size(verts, 1)
	ev = SparseArrays.spzeros(Int8, edge_num, verts_num)
	for i in 1:edge_num
		ev[i, alphas[alphas_keys[i]]] = 1
		ev[i, alphas[alphas_keys[i+1]]] = 1
	end
	return verts, ev
end


"""
	intersect_edges(
		V::Lar.Points,
		edge1::Lar.Cell,
		edge2::Lar.Cell
	)::Array{Tuple{Lar.Points, Float64}, 1}

Finds the intersection points (if there exist) between the two given edges.

Compute the new points over `edge1` given by the intersection with `edge2`.
For each intersection point it evaluates the scalar corresponding to the
relative distance from the first vertex of the edge.
Note that vertices of `edge1` are therefore never detected (see example 2).

See also: [`CAGD.frag_edge`](@ref)
---
# Examples
```jldoctest
# Cross
julia> V = [1.0 0.0; 0.0 1.0; 0.0 0.5; 0.5 1.0];
julia> copEV = Lar.coboundary_0([[1, 2], [3, 4]]);
julia> CAGD.intersect_edges(V, copEV[1, :], copEV[2, :])
1-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.25 0.75], 0.75)
```

```jldoctest
# Collinear
julia> V = [1.0 0.0; 0.0 1.0; 0.75 0.25; 0.5 0.5];
julia> copEV = Lar.coboundary_0([[1, 2], [3, 4]]);
julia> CAGD.intersect_edges(V, copEV[1, :], copEV[2, :])
2-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.75 0.25], 0.25)
 ([0.5 0.5], 0.5)
julia> CAGD.intersect_edges(V, copEV[2, :], copEV[1, :])
0-element Array{Tuple{Array{T,2} where T,Float64},1}
```
"""
function intersect_edges(
		V::Lar.Points,
		edge1::Lar.Cell,
		edge2::Lar.Cell
	)::Array{Tuple{Lar.Points, Float64}, 1}

	err = 10e-8

	x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
	x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
	ret = Array{Tuple{Lar.Points, Float64}, 1}()

	v1 = [x2-x1, y2-y1]
	v2 = [x4-x3, y4-y3]
	v3 = [x3-x1, y3-y1]
	ang1 = dot(LinearAlgebra.normalize(v1), LinearAlgebra.normalize(v2))
	ang2 = dot(LinearAlgebra.normalize(v1), LinearAlgebra.normalize(v3))
	parallel = 1-err < abs(ang1) < 1+err
	colinear = parallel && (1-err < abs(ang2) < 1+err || -err < norm(v3) < err)
	if colinear
		o = [x1 y1]
		v = [x2 y2] - o
		alpha = 1/dot(v,v')
		ps = [x3 y3; x4 y4]
		for i in 1:2
			a = alpha*dot(v',(reshape(ps[i, :], 1, 2)-o))
			if 0 < a < 1
				push!(ret, (ps[i:i, :], a))
			end
		end
	elseif !parallel
		denom = (v2[2])*(v1[1]) - (v2[1])*(v1[2])
		a = ((v2[1])*(-v3[2]) - (v2[2])*(-v3[1])) / denom
		b = ((v1[1])*(-v3[2]) - (v1[2])*(-v3[1])) / denom

		if -err < a < 1+err && -err <= b <= 1+err
			p = [(x1 + a*(x2-x1))  (y1 + a*(y2-y1))]
			push!(ret, (p, a))
		end
	end
	return ret
end


"""
merge_vertices!(
	model::CAGD.Model,
	[ edge_map::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}() ],
	[ err::Float64 = 1e-4 ]
)::Nothing

Modify `model` by compacting its points closer than `err` in a single one.

This method check one at time each vertex ``v`` in `model.G` and identifies
each other vertex within `err` with ``v`` itself.
The edges cochain `model.T[1]` is coherently modified (multiple edges between
two vertices are not allowed).
If an `edge_map` is given in input (this could be usefull during the planar
arrangement), then also the map is coherently modified and given back in output.

See also: [`CAGD.planar_arrangement_1`](@ref)
---
# Examples
```jldoctest
julia> model = CAGD.Model([
	0.5 0.0 0.5 1.0 0.5 1.0
	0.5 0.0 0.5 1.0 0.5 1.0
]);
julia> CAGD.addModelCells!(model, 1, [[1, 4], [3, 2], [5, 6], [1, 6], [5, 3]]);
julia> CAGD.merge_vertices!(model)

julia> model.G
2×3 Array{Float64,2}:
 0.5  0.0  1.0
 0.5  0.0  1.0

julia> Matrix(model.T[1])
2×3 Array{Int8,2}:
 1  0  1
 1  1  0
```
"""
function merge_vertices!(
		model::CAGD.Model,
		edge_map::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}(),
		err::Float64 = 1e-4
	)::Nothing
	V = convert(Lar.Points, model.G')
	EV = model.T[1]
	vertsnum = size(V, 1)
	edgenum = size(EV, 1)
	newverts = zeros(Int, vertsnum)
	# KDTree constructor needs an explicit array of Float64
	V = Array{Float64,2}(V)
	kdtree = KDTree(permutedims(V))

	# merge congruent vertices
	todelete = []
	i = 1
	for vi in 1:vertsnum
		if !(vi in todelete)
			nearvs = Lar.inrange(kdtree, V[vi, :], err)
			newverts[nearvs] .= i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i = i + 1
		end
	end
	nV = V[setdiff(collect(1:vertsnum), todelete), :]

	# merge congruent edges
	edges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
	oedges = Array{Tuple{Int, Int}, 1}(undef, edgenum)
	for ei in 1:edgenum
		v1, v2 = EV[ei, :].nzind
		edges[ei] = Tuple{Int, Int}(sort([newverts[v1], newverts[v2]]))
		oedges[ei] = Tuple{Int, Int}(sort([v1, v2]))
	end
	nedges = union(edges)
	nedges = filter(t->t[1]!=t[2], nedges)
	nedgenum = length(nedges)
	nEV = spzeros(Int8, nedgenum, size(nV, 1))
	# maps pairs of vertex indices to edge index
	etuple2idx = Dict{Tuple{Int, Int}, Int}()
	# builds `edge_map`
	for ei in 1:nedgenum
		nEV[ei, collect(nedges[ei])] .= 1
		etuple2idx[nedges[ei]] = ei
	end
	if !isempty(edge_map)
		for i in 1:length(edge_map)
			row = edge_map[i]
			row = map(x->edges[x], row)
			row = filter(t->t[1]!=t[2], row)
			row = map(x->etuple2idx[x], row)
			edge_map[i] = row
		end
	end
	# return new vertices and new edges
	model.G = convert(Lar.Points, nV')
	model.T[1] = nEV
	model.T[2] = SparseArrays.spzeros(Int8, 0, size(model, 1, 1))
	return
end


#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - CLEAN DECOMPOSITION
#-------------------------------------------------------------------------------

"""
	cleandecomposition(
		model::CAGD.Model,
		sigma::Union{Lar.Chain, Array{Int,1}},
		[ edge_map::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}() ]
	)

This function clears the `model` from all edges outside ``σ``.

This function take an arranged-edge `model` and a `Lar.Chain` ``σ`` made of
edges indices of the 1-chains `model.T[1]`.
The method drops from the Topology all the edges that are outside ``σ`` and
gives back the dropped model.
If an `edge_map` is given as an input it is also purged from the edges outside
w.r.t. ``σ`` and it is given as an output as well.

See also: [`CAGD.planar_arrangement`](@ref)
---
# Examples
```jldoctest
# Nested Triangles
julia> model = CAGD.Model([
	0.0  2.0  4.0  1.0  3.0  2.0  2.0 -2.0  6.0
	0.0  0.0  0.0  1.5  1.5  3.0 -3.0  3.0  3.0
]);
julia> CAGD.addModelCells!(model, 1, [
	[1, 3], [1, 6], [3, 6], [2, 4], [2, 5], [4, 5], [7, 8], [7, 9], [8, 9]
]);

julia> σ = SparseArrays.sparse([1; 1; 1; 0; 0; 0; 0; 0; 0]);

julia> model = CAGD.cleandecomposition(model, convert(Lar.Chain, σ))

julia> model.G
2×6 Array{Float64,2}:
 0.0  2.0  4.0  1.0  3.0  2.0
 0.0  0.0  0.0  1.5  1.5  3.0

julia> Matrix(model.T[1])
6×6 Array{Int8,2}:
 1  0  1  0  0  0
 1  0  0  0  0  1
 0  0  1  0  0  1
 0  1  0  1  0  0
 0  1  0  0  1  0
 0  0  0  1  1  0
```
"""
function cleandecomposition(
		model::CAGD.Model,
		sigma::Union{Lar.Chain, Array{Int,1}},
		edge_map::Array{Array{Int64,1},1} = Array{Array{Int64,1},1}()
	)::Union{CAGD.Model, Tuple{CAGD.Model, Array{Array{Int64,1},1}}}

	model = deepcopy(model)
	edge_map = deepcopy(edge_map)
	# Model point extraction to use Lar.point_in_face                           ##
	V = convert(Lar.Points, model.G')

	if issparse(sigma)
		sigma = sigma.nzind
	end
	todel = Array{Int,1}()
	sigma_edges = model.T[1][sigma, :]
	for e in 1 : model.T[1].m
		# Since the model is already arranged an edge can only be a sigma-edge,
		#  an intersection-free innner edge w.r.t. sigma face (not to purge),
		#  or an intersection-free outer edge (has to be deleted).
		if !(e in sigma)
			v1, v2 = CAGD.getModelCellVertices(model, 1, e)
			centroid = .5*(v1 + v2)

			if ! CAGD.point_in_face(centroid, V, sigma_edges)
				push!(todel, e)
			end
		end
	end

	CAGD.deleteModelCells!(model, 1, todel)
	CAGD.modelPurge!(model)

	if isempty(edge_map)
		return model
	end

	# Edges in edge_map must be updated too according to the same approach
	#  edges are deleted from models.
	for i in reverse(todel)
		for row in edge_map

			filter!(x->x!=i, row)

			for j in 1:length(row)
				if row[j] > i
					row[j] -= 1
				end
			end
		end
	end
	return model, edge_map
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - BICONNECTED COMPONENTS
#-------------------------------------------------------------------------------

"""
	biconnected_components(EV::Lar.ChainOp)

Compute sets of 2-cells on the same 3-cells biconnected components.

The method evaluates the 2-cells ``σ_i`` wich lies on the same 3-cells
biconnected component ``χ_j`` and gives back an array that contains an array
of ``σ_i`` for each ``j``.
Do note that if the cochain `EV` contains more copies of the same 2-cell then
it will be considered like a 3-cell.

Intermediate part of Planar Arrangement's algorithmic pipeline.

See also:
 - [`CAGD.planar_arrangement`](@ref) for the complete pipeline.
 - [`CAGD.planar_arrangement_2`](@ref).
---

# Examples
```jldoctest
julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0 0 0] #1 -> 1,2  |
		[1 0 1 0 0 0] #2 -> 1,3  |
		[1 0 0 1 0 0] #3 -> 1,4   |
		[1 0 0 0 1 0] #4 -> 1,5   |
		[1 0 0 0 0 1] #5 -> 1,6
		[0 1 1 0 0 0] #6 -> 2,3  |
		[0 0 0 1 1 0] #7 -> 4,5   |
	]));
julia> CAGD.biconnected_components(copEV)
2-element Array{Array{Int64,1},1}:
 [2, 6, 1]
 [4, 7, 3]
```

```jldoctest
julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0] #1 -> 1,2  |
		[1 1 0] #2 -> 1,2  |
		[1 0 1] #3 -> 1,2
	]));
julia> CAGD.biconnected_components(copEV)
1-element Array{Array{Int64,1},1}:
 [2, 1]
```
"""
function biconnected_components(EV::Lar.ChainOp)

	ps = Array{Tuple{Int, Int, Int}, 1}()
	es = Array{Tuple{Int, Int}, 1}()
	todel = Array{Int, 1}()
	visited = Array{Int, 1}()
	bicon_comps = Array{Array{Int, 1}, 1}()
	hivtx = 1

	function an_edge(point) # TODO: fix bug
		# error? : BoundsError: attempt to access 0×0 SparseMatrix ...
		edges = setdiff(EV[:, point].nzind, todel)
		if length(edges) == 0
			edges = [false]
		end
		edges[1]
	end

	function get_head(edge, tail)
		setdiff(EV[edge, :].nzind, [tail])[1]
	end

	function v_to_vi(v)
		i = findfirst(t->t[1]==v, ps)
		# seems findfirst changed from 0 to Nothing
		if typeof(i) == Nothing
			return false
		elseif i == 0
			return false
		else
			return ps[i][2]
		end
	end

	push!(ps, (1,1,1))
	push!(visited, 1)
	exit = false
	while !exit
		edge = an_edge(ps[end][1])
		if edge != false
			tail = ps[end][2]
			head = get_head(edge, ps[end][1])
			hi = v_to_vi(head)
			if hi == false
				hivtx += 1
				push!(ps, (head, hivtx, ps[end][2]))
				push!(visited, head)
			else
				if hi < ps[end][3]
					ps[end] = (ps[end][1], ps[end][2], hi)
				end
			end
			push!(es, (edge, tail))
			push!(todel, edge)
		else
			if length(ps) == 1
				found = false
				pop!(ps)
				for i in 1:size(EV,2)
					if !(i in visited)
						hivtx = 1
						push!(ps, (i, hivtx, 1))
						push!(visited, i)
						found = true
						break
					end
				end
				if !found
					exit = true
				end

			else
				if ps[end][3] == ps[end-1][2]
					edges = Array{Int, 1}()
					while true
						edge, tail = pop!(es)
						push!(edges, edge)
						if tail == ps[end][3]
							if length(edges) > 1
								push!(bicon_comps, edges)
							end
							break
						end
					end

				else
					if ps[end-1][3] > ps[end][3]
						ps[end-1] = (ps[end-1][1], ps[end-1][2], ps[end][3])
					end
				end
				pop!(ps)
			end
		end
	end
	bicon_comps = sort(bicon_comps, lt=(x,y)->length(x)>length(y))
	return bicon_comps
end


function monoconnected_components(EV::Lar.ChainOp)::Array{Array{Int,1},1}
	mono = Array{Array{Int,1},1}()
	tocheck = collect(1 : EV.m)
	while(!isempty(tocheck))
		vs = Array{Int,1}()
		vs_i = Array{Int,1}()
		vs_o = Array{Int,1}()
		evs = [tocheck[1]]
		evs_i = [tocheck[1]]
		evs_o = Array{Int,1}()

		while(!isempty(vs_i) || !isempty(evs_i))
			if !isempty(evs_i)
				e = evs_i[1]
				push!(evs_o, e);
				push!(vs, EV[e, :].nzind...)
			end
			vs_i = setdiff(vs, vs_o)
			if !isempty(vs_i)
				v = vs_i[1]
				push!(vs_o, v);
				push!(evs, EV[:, v].nzind...)
			end
			evs_i = setdiff(evs, evs_o)
		end

		push!(mono, unique(evs))
	end
	return mono
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - REMOVE MAPPING DANGLING EDGES
#-------------------------------------------------------------------------------

"""
	remove_mapping_dangling_edges(
		edge_map::Array{Array{Int64,1},1},
		bicon_comps::Array{Array{Int64,1},1}
	)::Array{Array{Int64,1},1}

Remove from the edge_map the edges not belonging to a biconnected component.

This utility function deletes from the `edge_map` all the edges that do not
belongs to a biconnected component by ordinatelly recompatting the indices of
the remaning edges.

See also: [`CAGD.planar_arrangement`](@ref).
"""
function remove_mapping_dangling_edges(
		edge_map::Array{Array{Int64,1},1},
		bicon_comps::Array{Array{Int64,1},1}
	)::Array{Array{Int64,1},1}

	edge_map = deepcopy(edge_map)

	edges = sort(union(bicon_comps...))
								  #size(model.T[1, 1]) == max(edge_map...)
	todel = sort(setdiff(collect(1:max(max(edge_map...)...)), edges))

	for i in reverse(todel)
		for row in edge_map

			filter!(x->x!=i, row)

			for j in 1:length(row)
				if row[j] > i
					row[j] -= 1
				end
			end
		end
	end

	return edge_map
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 METHODS
#    componentgraph (aka TGW) from tgw2.jl
#    cell_merging
#-------------------------------------------------------------------------------

"""
	cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

Cells composing for the Topological Gift Wrapping algorithm.

This is the online part of the TGW algorithm.

See also: [`CAGD.planar_arrangement_2`](@ref).
"""
function cell_merging(
		n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
	)
	function bboxes(V::Lar.Points, indexes::Lar.ChainOp)
		boxes = Array{Tuple{Any, Any}}(undef, indexes.n)
		for i in 1:indexes.n
			v_inds = indexes[:, i].nzind
			boxes[i] = Lar.bbox(V[v_inds, :])
		end
		boxes
	end
	# initialization
	sums = Array{Tuple{Int, Int, Int}}(undef, 0)
	# assembling child components with father components
	for father in 1:n
		if sum(containment_graph[:, father]) > 0
			father_bboxes = bboxes(
				V, abs.(EVs[father]')*abs.(boundaries[father]')
			)
			for child in 1:n
				if containment_graph[child, father] > 0
					child_bbox = shell_bboxes[child]
					for b in 1:length(father_bboxes)
						if Lar.bbox_contains(father_bboxes[b], child_bbox)
							push!(sums, (father, b, child))
							break
						end
					end
				end
			end
		end
	end
	# offset assembly initialization
	EV = vcat(EVs...)
	edgenum = size(EV, 1)
	facenum = sum(map(x->size(x,1), boundaries))
	FE = spzeros(Int8, facenum, edgenum)
	shells2 = spzeros(Int8, length(shells), edgenum)
	r_offsets = [1]
	c_offset = 1
	# submatrices construction
	for i in 1:n
		min_row = r_offsets[end]
		max_row = r_offsets[end] + size(boundaries[i], 1) - 1
		min_col = c_offset
		max_col = c_offset + size(boundaries[i], 2) - 1
		FE[min_row:max_row, min_col:max_col] = boundaries[i]
		shells2[i, min_col:max_col] = shells[i]
		push!(r_offsets, max_row + 1)
		c_offset = max_col + 1
	end
	# offsetting assembly of component submatrices
	for (f, r, c) in sums
		FE[r_offsets[f]+r-1, :] += shells2[c, :]
	end

	return EV, FE
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1
#-------------------------------------------------------------------------------


"""
	planar_arrangement_1(
		model::CAGD.Model,
		[sigma::Lar.Chain = spzeros(Int8, 0)],
		[return_edge_map::Bool = false]
	)

First part of arrangement's algorithmic pipeline.

This function computes the pairwise intersection between each edge of a given
2D cellular complex 1-skeleton. The computation is speeded up via the evaluation
of the Spatial Index. See [`CAGD.spaceindex`](@ref).

See also: [`CAGD.planar_arrangement`](@ref) for the complete pipeline.
It uses:
 - [`CAGD.frag_edges_channel`](@ref)
 - [`CAGD.frag_edges`](@ref)
 - [`CAGD.merge_vertices!`](@ref)
---
# Examples
```jldoctest
julia> model = CAGD.Model([
	0.0 0.5 0.0 0.5 0.3 1.0 0.3 1.0;
	0.0 0.0 1.0 1.0 0.5 0.5 1.0 1.0
]);
julia> CAGD.addModelCells!(model, 1, [
	[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [5, 7], [6, 8]
]);
julia> (model, edge_map) =
	CAGD.planar_arrangement_1(model, spzeros(Int8, 0), true);
julia> model.G
2×9 Array{Float64,2}:
 0.0  0.5  0.0  0.5  0.3  0.5  0.3  1.0  1.0
 0.0  0.0  1.0  1.0  1.0  0.5  0.5  0.5  1.0

julia> Matrix(model.T[1])
11×9 Array{Int8,2}:
 1  1  0  0  0  0  0  0  0
 0  0  1  1  0  0  0  0  0
 0  0  0  1  1  0  0  0  0
 1  0  1  0  0  0  0  0  0
 0  1  0  0  0  1  0  0  0
 0  0  0  0  1  1  0  0  0
 0  0  0  0  0  1  1  0  0
 0  0  0  0  0  1  0  1  0
 0  0  0  0  1  0  0  0  1
 0  0  0  1  0  0  1  0  0
 0  0  0  0  0  0  0  1  1

julia> edge_map
8-element Array{Array{Int64,1},1}:
 [1]
 [2, 3]
 [4]
 [5, 6]
 [7, 8]
 [3, 9]
 [10]
 [11]

```
"""
function planar_arrangement_1(
		model::CAGD.Model,
		sigma::Lar.Chain = spzeros(Int8, 0),
		return_edge_map::Bool = false
	)::Union{
		CAGD.Model,
		Tuple{CAGD.Model, Array{Array{Int64,1},1}},
		Tuple{CAGD.Model, Array{Int,1}},
		Tuple{CAGD.Model, Array{Int,1}, Array{Array{Int64,1},1}}
	}
	# data structures initialization
	# V = convert(Lar.Points, model.G')
	# copEV = model.T[1]
	edgenum = size(model, 1, 1)
	edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
	rV = Lar.Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# spaceindex computation
	bigPI = CAGD.spaceIndex(model, 1)


	# sequential (iterative) processing of edge fragmentation
	for i in 1:edgenum
		v, ev = CAGD.frag_edge(model, i, bigPI)
		newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
		edge_map[i] = newedges_nums
		finalcells_num += size(ev, 1)
		rV = convert(Lar.Points, rV)
		rV, rEV = Lar.skel_merge(rV, rEV, v, ev)
	end

	# merging of close vertices and edges (2D congruence)
	model = CAGD.Model(convert(Lar.Points, rV'))
	CAGD.addModelCells!(model, 1, rEV)
	CAGD.merge_vertices!(model, edge_map)

	if isempty(sigma)
		if return_edge_map
			return model, edge_map
		else
			return model
		end
	end

	# Sigma update by mapping each of its edges into newly born edges
	new_sigma = Array{Int,1}()
	map(i->new_sigma=union(new_sigma, edge_map[i]), sigma.nzind)
	if return_edge_map
		return model, new_sigma, edge_map
	else
		return model, new_sigma
	end
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2
#-------------------------------------------------------------------------------

"""
	planar_arrangement_2(
		model::CAGD.Model,
		bicon_comps::{Array{Array,1},1}
	)::CAGD.Model

Second part of arrangement's algorithmic pipeline.

This function is the complete Topological Gift Wrapping (TGW) algorithm that
is firstly locally used in order to decompose the 2-cells and then globally
to generate the 3-cells of the arrangement of the ambient space ``E^3``.

During this process each dangling 1-cell is removed.
Do note that the isolated 2-cells are not removed by this procedure.

See also: [`CAGD.planar_arrangement`](@ref) for the complete pipeline.
It uses:
 - [`CAGD.biconnected_components`](@ref)
 - [`CAGD.componentgraph`](@ref)
 - [`CAGD.cell_merging`](@ref)
---

# Examples
```jldoctest
# Triforce
julia> model = CAGD.Model([
	0.0 0 2 0 0 1 4 4 3 2 4 4
	0.0 2 0 3 4 4 2 0 0 4 4 2
]);

julia> CAGD.addModelCells!(model, 1, [
	[ 1,  2], [ 2,  3], [ 3,  1],
	[ 4,  5], [ 5,  6], [ 6,  7], [ 7,  8], [ 8,  9], [ 9, 4],
	[10, 11], [11, 12], [12, 10]
])

julia> bicon_comps = CAGD.biconnected_components(model.T[1])
3-element Array{Array{Int64,1},1}:
 [9, 8, 7, 6, 5, 4]
 [3, 2, 1]
 [12, 11, 10]

julia> model = CAGD.planar_arrangement_2(model, bicon_comps);

julia> Matrix(model.T[2])
3×12 Array{Int8,2}:
 -1  -1  -1  -1  -1  1   0   0  0   0   0  0
  0   0   0   0   0  0  -1  -1  1   0   0  0
  0   0   0   0   0  0   0   0  0  -1  -1  1
```
"""
function planar_arrangement_2(
		model::CAGD.Model,
		bicon_comps::Array{Array{Int64,1},1}
	)::CAGD.Model

	#==
	test to remove
	@assert isequal(
		bicon_comps,
		CAGD.biconnected_components(model.T[1])
	)
	==#

	# component graph
	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes =
		CAGD.componentgraph(
			convert(Lar.Points, model.G'), model.T[1], bicon_comps
		)

	copEV, FE = CAGD.cell_merging(
	   	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
	)

	model = CAGD.Model(convert(Lar.Points, V'))
	CAGD.addModelCells!(model, 1, copEV)
	CAGD.addModelCells!(model, 2, FE)
	CAGD.modelPurge!(model)
	return model
end


#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE
#-------------------------------------------------------------------------------

"""
	planar_arrangement(
		model::CAGD.Model,
		[sigma::Lar.Chain],
		[return_edge_map::Bool]
	)::Union{CAGD.Model, Tuple{CAGD.Model, Array{Array{Int, 1}, 1}}}

Compute the arrangement on the given cellular complex 1-skeleton in 2D.

A cellular complex is arranged when the intersection of every possible pair
of cell of the complex is empty and the union of all the cells is the
whole Euclidean space.
The methods default input/output is a `CAGD.Model` with at least Geometry and
1-cell Topoly.
Mutiprocessing is by default disabled.

It uses:
 - [`CAGD.planar_arrangement_1`](@ref)
 - [`CAGD.cleandecomposition`](@ref)
 - [`CAGD.biconnected_components`](@ref)
 - [`CAGD.remove_mapping_dangling_edges`](@ref)
 - [`CAGD.planar_arrangement_2`](@ref)

## Additional arguments:
 - `sigma::Chain`: if specified, the method return only cells contained inside
		this cell. _Default_: empty cell.
 - `return_edge_map::Bool`: builds an `edge_map` which maps the input edges
		to the output's. _Defaults_: `false`.
"""
function planar_arrangement(
		model::CAGD.Model,
		sigma::Lar.Chain = spzeros(Int8, 0),                                    ## TEST
		return_edge_map::Bool = false
	)::Union{CAGD.Model, Tuple{CAGD.Model, Array{Array{Int, 1}, 1}}}

	# Check dimensional condition
	@assert length(model) == 2

	#planar_arrangement_1
	if isempty(sigma)
		model, edge_map =
			CAGD.planar_arrangement_1(model, sigma, true)
	else
		model, sigma, edge_map =
			CAGD.planar_arrangement_1(model, sigma, true)
	end

#   print("PLANAR_ARRANGEMENT_1 COMPLETE\n")

	# cleandecomposition
	if !isempty(sigma)
		model, edge_map = CAGD.cleandecomposition(model, sigma, edge_map)
#       print("DECOMPOSITION CLEANED\n")
	end

	# generation of biconnected components
	bicon_comps = CAGD.biconnected_components(model.T[1])
	#mocon_comps = CAGD.monoconnected_components(model.T[1])

#	print("BICONNECTED COMPONENTS FOUND\n")

	if isempty(bicon_comps)
		println("No biconnected components found.\n")
		if (return_edge_map)
			return (model(), nothing)
		else
			return model()
		end
	end

	# Removing Dangling edges from edge map
	if return_edge_map
		@assert size(model, 1, 1) == max((edge_map...)...)                           # CTR
		edge_map = CAGD.remove_mapping_dangling_edges(edge_map, bicon_comps)
	end

#	print("DANGLING EDGES REMOVED\n")

	# Planar_arrangement_2
	model = CAGD.planar_arrangement_2(model, bicon_comps)

#	print("PLANAR ARRANGEMENT 2 COMPLETE\n")

	if size(model, 0, 2) - size(model, 1, 1) + size(model, 2, 1) != 1
		@show size(model, 0, 2) - size(model, 1, 1) + size(model, 2, 1)
	end

	if return_edge_map
		 return model, edge_map
	else
		 return model
	end
end


"""
	point_in_face(point, V::Points, copEV::ChainOp)
Check if `point` is not outside the area of the face bounded by the edges in `copEV`
"""
function point_in_face(point, V::Lar.Points, copEV::Lar.ChainOp)

	function pointInPolygonClassification(V,EV)

		function crossingTest(new, old, status, count)
			if status == 0
				status = new
				return status, (count + 0.5)
			else
				if status == old
					return 0, (count + 0.5)
				else
					return 0, (count - 0.5)
				end
			end
		end

		# Outputs a tileCoder that takes as input a point and classify its
		#  octant position with respect to the supercrescent list N-S-E-W:
		#  1001 | 0001 | 0101
		#  -----+------+-----
		#  1000 | 0000 | 0100
		#  -----+------+-----
		#  0101 | 0010 | 0110
		function setTile(box)
			tiles = [[9,1,5],[8,0,4],[10,2,6]]
			b1,b2,b3,b4 = box
			function tileCode(point)
				x,y = point
				code = 0
				if y>b1 code=code|1 end
				if y<b2 code=code|2 end
				if x>b3 code=code|4 end
				if x<b4 code=code|8 end
				return code
			end
			return tileCode
		end

		function pointInPolygonClassification0(pnt)
			# Set the tile coder w.r.t. the point we are checking
			x,y = pnt
			xmin,xmax,ymin,ymax = x,x,y,y
			tilecode = setTile([ymax,ymin,xmax,xmin])
			count,status = 0,0

			for k in 1:EV.m
				edge = EV[k,:]
				p1, p2 = V[edge.nzind[1], :], V[edge.nzind[2], :]
				(x1,y1),(x2,y2) = p1,p2
				c1,c2 = tilecode(p1),tilecode(p2)
				# For each edge we identify via:
				#  - xor: the different sides
				#  -  or: all the direction involved
				#  - and: the cardinal direction (if unique)
				c_edge, c_un, c_int = xor(c1, c2), c1|c2, c1&c2

				if (c_edge == 0) & (c_un == 0) return "p_on"
				elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
				elseif c_edge == 3
					if c_int == 0 return "p_on"
					elseif c_int == 4 count += 1 end
				elseif c_edge == 15
					x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2
					if x_int > x count += 1
					elseif x_int == x return "p_on" end
				elseif (c_edge == 13) & ((c1==4) | (c2==4))
						status, count = crossingTest(1,2,status,count)
				elseif (c_edge == 14) & ((c1==4) | (c2==4))
						status, count = crossingTest(2,1,status,count)
				elseif c_edge == 7 count += 1
				elseif c_edge == 11 count = count
				elseif c_edge == 1
					if c_int == 0 return "p_on"
					elseif c_int == 4
						status, count = crossingTest(1,2,status,count)
					end
				elseif c_edge == 2
					if c_int == 0 return "p_on"
					elseif c_int == 4
						status, count = crossingTest(2,1,status,count)
					end
				elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
				elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
				elseif c_edge == 5
					if (c1==0) | (c2==0) return "p_on"
					else
						status, count = crossingTest(1,2,status,count)
					end
				elseif c_edge == 6
					if (c1==0) | (c2==0) return "p_on"
					else
						status, count = crossingTest(2,1,status,count)
					end
				elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
				elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
				end
			end

			if (round(count)%2)==1
				return "p_in"
			else
				return "p_out"
			end
		end
		return pointInPolygonClassification0
	end

	ret = pointInPolygonClassification(V, copEV)(point)
	return "p_in" == ret || ret == "p_on"
end
