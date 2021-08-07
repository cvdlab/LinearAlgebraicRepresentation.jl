"""
    componentgraph(V, copEV, bicon_comps)

Topological Gift Wrapping algorithm on 2D skeletons.

This is the offline part of the TGW algorithm. It takes in input a model and its
biconnected components mapping and evaluates usefull informations:
 1. Number of biconnected components.
 2. Component Graph of the biconnected structure.
 3. The 1-cells structure (UNMODIFIED).                                         # Could be removed?
 4. Association between non-dangling 2-cells and their orientation (∀ component)
 5. Association between 3-cells and 2-cells (with orientation, ∀ component).
 6. Association between 3-cells and their orientation (∀ component).
 7. Shell bounding boxes of the components.

See also: [`CAGD.planar_arrangement_2`](@ref).
It uses:
 - [`CAGD.get_external_cycle`](@ref)
 - [`CAGD.pre_containment_test`](@ref)
 - [`CAGD.prune_containment_graph`](@ref)
 - [`CAGD.transitive_reduction!`](@ref)
---
# Examples
```jldoctest
# Papillon
julia> V = [0.0 0.0; 0.0 3.0; 2.0 1.5; 4.0 0.0; 4.0 3.0];

julia> copEV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0 0] #1 -> 1,2
		[0 1 1 0 0] #2 -> 2,3
		[1 0 1 0 0] #3 -> 3,1
		[0 0 1 1 0] #4 -> 3,4
		[0 0 0 1 1] #5 -> 4,5
		[0 0 1 0 1] #6 -> 3,5
	]));

julia> copFE = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 1 0 0 0] #1 -> 1,2,3
		[0 0 0 1 1 1] #2 -> 4,5,6
		[1 1 1 1 1 1] #3 -> 1,2,3,4,5,6 External
	]));

julia> bicon_comps = [[1, 2, 3], [4, 5, 6]];

julia> CAGD.componentgraph(V, copEV, bicon_comps)[1]
2

julia> CAGD.componentgraph(V, copEV, bicon_comps)[2]
2×2 SparseMatrixCSC{Int8,Int64} with 0 stored entries

julia> CAGD.componentgraph(V, copEV, bicon_comps)[4]
2-element Array{SparseMatrixCSC{Int8,Int64},1}:

  [1, 1]  =  -1
  [3, 1]  =  -1
  [1, 2]  =  1
  [2, 2]  =  -1
  [2, 3]  =  1
  [3, 3]  =  1

  [1, 3]  =  -1
  [3, 3]  =  -1
  [1, 4]  =  1
  [2, 4]  =  -1
  [2, 5]  =  1
  [3, 5]  =  1

julia> CAGD.componentgraph(V, copEV, bicon_comps)[5]
2-element Array{SparseMatrixCSC{Int8,Int64},1}:

  [1, 1]  =  -1
  [1, 2]  =  -1
  [1, 3]  =  1

  [1, 1]  =  1
  [1, 2]  =  1
  [1, 3]  =  -1

julia> CAGD.componentgraph(V, copEV, bicon_comps)[6]
2-element Array{SparseVector{Int8,Int64},1}:
   [1]  =  1
  [2]  =  1
  [3]  =  -1
   [1]  =  -1
  [2]  =  -1
  [3]  =  1

julia> CAGD.componentgraph(V, copEV, bicon_comps)[7]
2-element Array{Any,1}:
 ([0.0 0.0], [2.0 3.0])
 ([2.0 0.0], [4.0 3.0])
```
"""
function componentgraph(V, copEV, bicon_comps)
    # arrangement of isolated components
	n = size(bicon_comps, 1)
   	shells = Array{Lar.Chain, 1}(undef, n)
	boundaries = Array{Lar.ChainOp, 1}(undef, n)
	EVs = Array{Lar.ChainOp, 1}(undef, n)
    # for each component
	for p in 1 : n
		ev = copEV[sort(bicon_comps[p]), :]
        # computation of 2-cells
		fe = Lar.Arrangement.minimal_2cycles(V, ev)
        # exterior cycle
		shell_num = CAGD.get_external_cycle(V, ev, fe)
        # decompose each fe (co-boundary local to component)
		EVs[p] = ev
		tokeep = setdiff(1:fe.m, shell_num)
		boundaries[p] = fe[tokeep, :]
		shells[p] = fe[shell_num, :]
    end
    # computation of bounding boxes of isolated components
	shell_bboxes = []
	for i in 1 : n
    	vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
   		push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
	end
    # computation and reduction of containment graph
	containment_graph = CAGD.pre_containment_test(shell_bboxes)
	containment_graph = CAGD.prune_containment_graph(V,EVs,shells,containment_graph)
	CAGD.transitive_reduction!(containment_graph)
	return n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
end


#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 - COMPONENT GRAPH (TGW) METHODS
#    get_external_cycle
#    pre_containment_test
#    prune_containment_graph
#    transitive_reduction!
#-------------------------------------------------------------------------------

"""
    get_external_cycle(model::CAGD.Model)::Union{Int64, Nothing}
    get_external_cycle(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)

Evaluates the index of the external 3D-cell of a Topological Model.

This method looks for and retrieve the external cycle of the `model`.
If the cell does not exist then it returns `nothing`.
The method is also callable by the signature `(V, EV, FE)`.

See also: [`CAGD.componentgraph`](@ref).
---
# Examples
```jldoctest
# ``K^4`` graph
julia> model = CAGD.Model(hcat([[0.0,0.0], [4.0,0.0], [2.0,3.0], [2.0,1.5]]...));
julia> CAGD.addModelCells!(model, 1, [[1,2],[2,3],[3,1],[1,4],[2,4],[3,4]]);
julia> CAGD.addModelCells!(model, 2, [[1,4,5],[2,5,6],[3,4,6],[1,2,3]]);

julia> CAGD.get_external_cycle(model)
4

julia> CAGD.deleteModelCell!(model, 2, 4);
julia> typeof(CAGD.get_external_cycle(model))
Nothing
```
"""
function get_external_cycle(model::CAGD.Model)::Union{Int64, Nothing}
    CAGD.get_external_cycle(
        convert(Lar.Points, model.G'),
        model.T[1],
        model.T[2]
    )
end

function get_external_cycle(
        V::Lar.Points,
        EV::Lar.ChainOp,
        FE::Lar.ChainOp
    )::Union{Int64, Nothing}

    FV = abs.(FE)*EV

    vs = sparsevec(mapslices(sum, abs.(EV), dims=1)').nzind
    minv_x1 = maxv_x1 = minv_x2 = maxv_x2 = pop!(vs)
    for i in vs
        if V[i, 1] > V[maxv_x1, 1]
            maxv_x1 = i
        elseif V[i, 1] < V[minv_x1, 1]
            minv_x1 = i
        end
        if V[i, 2] > V[maxv_x2, 2]
            maxv_x2 = i
        elseif V[i, 2] < V[minv_x2, 2]
            minv_x2 = i
        end
    end
    cells = intersect(
        FV[:, minv_x1].nzind,
        FV[:, maxv_x1].nzind,
        FV[:, minv_x2].nzind,
        FV[:, maxv_x2].nzind
    )
    if length(cells) == 1
        return cells[1]
    else
        for c in cells
            if Lar.face_area(V, EV, FE[c, :]) < 0
                return c
            end
        end
    end
end

"""
    pre_containment_test(
        bboxes::Array{Tuple{Array{Float64,2},Array{Float64,2}},1}
    )::SparseMatrixCSC{Int8,Int64}

Generate the containment graph associated to `bboxes`.

The function evaluate a `SparseArrays{Int8, n, n}` PRE-containmnet graph meaning
it stores `1` when the bbox associated to the `i`-th component contains entirely
the bbox associated to the `j`-th one. (where `n = length(bboxes)`)

See also: [`CAGD.componentgraph`](@ref).

---

# Examples
```jldoctest
# Planar Tiles
julia> bboxes = [
		([0.0 0.0], [1.0 1.0])
		([0.0 0.0], [0.4 0.4])
		([0.6 0.0], [1.0 0.4])
		([0.0 0.6], [0.4 1.0])
		([0.6 0.6], [1.0 1.0])
		([2.0 3.0], [2.0 3.0]) # no intersection
	];

julia> CAGD.pre_containment_test(bboxes)
6x6 SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [2, 1]  =  1
  [3, 1]  =  1
  [4, 1]  =  1
  [5, 1]  =  1
```
"""
function pre_containment_test(
        bboxes#::Array{Tuple{Array{Float64,2},Array{Float64,2}},1}              #
    )::SparseMatrixCSC{Int8,Int64}
    n = length(bboxes)
    containment_graph = spzeros(Int8, n, n)

    for i in 1:n
        for j in 1:n
            if i != j && Lar.bbox_contains(bboxes[j], bboxes[i])
                containment_graph[i, j] = 1
            end
        end
    end

    return containment_graph
end

"""
    prune_containment_graph(
            V::Lar.Points,
            EVs::Array{Lar.ChainOp, 1},
            shells::Array{Lar.Chain,1},
            graph::SparseMatrixCSC{Int8, Int64}
        )::SparseMatrixCSC{Int8, Int64}

Prunes the containment `graph` from the non-included faces.

This method prunes the containment graph eliminating the included bouning boxes
that do not corresponds to included faces.

Do note that this method expects to work on biconnectet clusters of faces.

See also: [`CAGD.componentgraph`](@ref).
```
"""
function prune_containment_graph(
        V::Lar.Points,
        EVs::Array{Lar.ChainOp, 1},
        shells::Array{Lar.Chain,1},
        graph::SparseMatrixCSC{Int8, Int64}
    )::SparseMatrixCSC{Int8, Int64}

    graph = deepcopy(graph) # copy needed
    n = length(EVs)
    length(shells) == n || throw(ArgumentError("EVs and shells not coherent"))
    size(graph) == (n, n) || throw(ArgumentError("graph dim not coherent"))

    for i in 1 : n
        # Take a point of the i-th shell by taking a vertex of its first edge
        an_edge = shells[i].nzind[1]
        origin_index = EVs[i][an_edge, :].nzind[1]
        origin = V[origin_index, :]

        # Check the containment for every other cell by testing if such a point
        #   is inner to the face defined byt that cell
        for j in 1 : n
            if i != j && graph[i, j] == 1
                shell_edge_indexes = shells[j].nzind
                ev = EVs[j][shell_edge_indexes, :]

                if !CAGD.point_in_face(origin, V, ev)
                    graph[i, j] = 0
                end
            end
         end

     end
     return dropzeros(graph)
end

"""
    transitive_reduction!(graph)

Prune the redundancies in containment `graph`.

The ``m``-layer containment `graph` is reduced to the corresponding 1-layer
containment `graph` where only the first level nested components are considered.

See also: [`CAGD.componentgraph`](@ref).

---

# Examples
```jldoctest
julia> containment = [
		0 0 0 0 0 0
		0 0 0 0 0 0
		1 0 0 0 0 0
		1 0 0 0 0 0
		1 0 1 0 0 0
		1 0 1 0 1 0
	];

julia> CAGD.transitive_reduction!(containment)

julia> containment
6×6 Array{Int64,2}:
 0  0  0  0  0  0
 0  0  0  0  0  0
 1  0  0  0  0  0
 1  0  0  0  0  0
 0  0  1  0  0  0
 0  0  0  0  1  0
```
"""
function transitive_reduction!(graph)
    n = size(graph, 1)
    for j in 1:n
        for i in 1:n
            if graph[i, j] > 0
                for k in 1:n
                    if graph[j, k] > 0
                        graph[i, k] = 0
                    end
                end
            end
        end
    end
end
