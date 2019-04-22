using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation


#--------------------------------------------------------------------------------------------------------------------------------
#                               PIPELINE PART 1
#--------------------------------------------------------------------------------------------------------------------------------

"""
    frag_edge_channel(in_chan, out_chan, V, EV, bigPI)

Utility function for parallel edge fragmentation.

This function handles the edge fragmentation in the first part of arrangement's algorithmic pipeline
(see also [`Lar.Arrangement.planar_arrangement_1`](@ref)) during multiprocessing computation.
In order to do so it needs two `Distributed.RemoteChannel`s, one with the inputs and one for outputs.

See also: [`Lar.Arrangement.frag_edge`](@ref)
"""
function frag_edge_channel( in_chan::Distributed.RemoteChannel{Channel{Int64}},
                            out_chan::Distributed.RemoteChannel{Channel{Tuple}},
                            V::Lar.Points,
                            EV::Lar.ChainOp,
                            bigPI::Array{Array{Int64,1},1})
    run_loop = true
    while run_loop
        edgenum = take!(in_chan)
        if edgenum != -1
            put!(out_chan, (edgenum, Lar.Arrangement.frag_edge(V, EV, edgenum, bigPI)))
        else
            run_loop = false
        end
    end
end


"""
	frag_edge(V, EV, edge_idx, bigPI)

Splits the `edge_idx`-th edge in `EV`.

This method splits the `edge_idx`-th edge in `EV` into several parts by confronting it
with the others that intersect its bounding box (see also [`Lar.Arrangement.intersect_edges`](@ref)).

The method returns a set of the new vertices that the segment is made of (with redundancies) and
the associated cochain (with no redundancies).

See also:
 - [`Lar.Arrangement.planar_arrangement_1`](@ref)
 - [`Lar.Arrangement.frag_edge_channel`](@ref)

---

# Examples
```jldoctest
julia> V = [1.0 0.0; 0.0 1.0; 0.0 0.5; 0.5 1.0; 1.0 1.0]; # By Rows!

julia> EV = [[1, 2], [2, 5], [3, 4], [4, 5]];

julia> cop_EV = Lar.coboundary_0(EV::Lar.Cells);

julia> bigPI = Lar.spaceindex((convert(Lar.Points, V'), EV));

julia> Lar.Arrangement.frag_edge(V, cop_EV, 1, bigPI)[1]
5×2 Array{Float64,2}:
 1.0   0.0 
 0.0   1.0 
 1.0   0.0 
 0.25  0.75
 0.0   1.0 

julia> Lar.Arrangement.frag_edge(V, cop_EV, 1, bigPI)[2]
2×5 SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [1, 1]  =  1
  [2, 2]  =  1
  [1, 4]  =  1
  [2, 4]  =  1
```
"""
function frag_edge(V::Lar.Points, EV::Lar.ChainOp, edge_idx::Int, bigPI::Array{Array{Int64,1},1})::Tuple{Lar.Points, Lar.ChainOp}
    alphas = Dict{Float64, Int}()
    edge = EV[edge_idx, :]
    verts = V[edge.nzind, :]
    for i in bigPI[edge_idx]
        if i != edge_idx && edge_idx in bigPI[i]    #<------------------------- Could be usefull?
            intersection = Lar.Arrangement.intersect_edges(V, edge, EV[i, :])
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
	intersect_edges(V, edge1, edge2)

Finds the intersection point (if there exists) between the two given edges.

This method compute the points where `edge2` intersect `edge1`.
If they are collinear only the vertices of `edge2` are considered (see the second example).

See also: [`Lar.Arrangement.frag_edge`](@ref)

---

# Examples
```jldoctest
julia> V = [1.0 0.0; 0.0 1.0; 0.0 0.5; 0.5 1.0]; # By Rows!

julia> EV = [[1, 2], [3, 4]];

julia> cop_EV = Lar.coboundary_0(EV::Lar.Cells);

julia> Lar.Arrangement.intersect_edges(V, cop_EV[1, :], cop_EV[2, :])
1-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.25 0.75], 0.75)
```

```jldoctest
julia> V = [1.0 0.0; 0.0 1.0; 0.75 0.25; 0.5 0.5]; # By Rows!

julia> EV = [[1, 2], [3, 4]];

julia> cop_EV = Lar.coboundary_0(EV::Lar.Cells);

julia> Lar.Arrangement.intersect_edges(V, cop_EV[1, :], cop_EV[2, :])
2-element Array{Tuple{Array{T,2} where T,Float64},1}:
 ([0.75 0.25], 0.25)
 ([0.5 0.5], 0.5) 

julia> Lar.Arrangement.intersect_edges(V, cop_EV[2, :], cop_EV[1, :])
0-element Array{Tuple{Array{T,2} where T,Float64},1}
```
"""
function intersect_edges(V::Lar.Points, edge1::Lar.Cell, edge2::Lar.Cell)::Array{Tuple{Lar.Points, Float64}, 1}
    err = 10e-8

    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Lar.Points, Float64}, 1}()

    v1 = [x2-x1, y2-y1];
    v2 = [x4-x3, y4-y3];
    v3 = [x3-x1, y3-y1];

    ang1 = Lar.dot(Lar.normalize(v1), Lar.normalize(v2))
    ang2 = Lar.dot(Lar.normalize(v1), Lar.normalize(v3))
    
    parallel = 1-err < abs(ang1) < 1+err
    colinear = parallel && (1-err < abs(ang2) < 1+err || -err < Lar.norm(v3) < err)
    

    if colinear
        o = [x1 y1] 
        v = [x2 y2] - o
        alpha = 1/Lar.dot(v,v')
        ps = [x3 y3; x4 y4]
        for i in 1:2
            a = alpha*Lar.dot(v',(reshape(ps[i, :], 1, 2)-o))
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
    merge_vertices!(V, EV[, edge_map[, err]])

Compact the vertices closer than `err` in a single one.

This method check one at time each vertex ``v`` in `V` and identifies each other vertex within `err` with ``v`` itself.
The cochain `EV` is coherently modified (multiple edges between two vertices are not allowed).
If an `edge_map` is given in input (this could be usefull during the planar arrangements), then also
the map is coherently modified and given back in output.

See also: [`Lar.Arrangement.planar_arrangement_1`](@ref)

---

# Examples
```jldoctest
julia> V = [0.5 0.5; 0.0 0.0; 0.5 0.5; 1.0 1.0; 0.5 0.5; 1.0 1.0]; # By Rows!

julia> EV = [[1, 4], [3, 2], [5, 6], [1, 6], [5, 3]];

julia> cop_EV = Lar.coboundary_0(EV::Lar.Cells);

julia> Lar.Arrangement.merge_vertices!(V, cop_EV)[1]
3×2 Array{Float64,2}:
 0.5  0.5
 0.0  0.0
 1.0  1.0

julia> Lar.Arrangement.merge_vertices!(V, cop_EV)[2]
2×3 SparseArrays.SparseMatrixCSC{Int8,Int64} with 4 stored entries:
  [1, 1]  =  1
  [2, 1]  =  1
  [2, 2]  =  1
  [1, 3]  =  1
```
"""
function merge_vertices!(V::Lar.Points, EV::Lar.ChainOp, edge_map::Array{Array{Int64,1},1}=[[-1]], err=1e-4)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    newverts = zeros(Int, vertsnum)
    # KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    kdtree = KDTree(permutedims(V))

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
    
    etuple2idx = Dict{Tuple{Int, Int}, Int}()
    
    for ei in 1:nedgenum
        nEV[ei, collect(nedges[ei])] .= 1
        etuple2idx[nedges[ei]] = ei
    end
    
    if edge_map != [[-1]]
        for i in 1:length(edge_map)
            row = edge_map[i]
            row = map(x->edges[x], row)
            row = filter(t->t[1]!=t[2], row)
            row = map(x->etuple2idx[x], row)
            edge_map[i] = row
        end
        return Lar.Points(nV), nEV, edge_map
    end
    

    return Lar.Points(nV), nEV
end



#--------------------------------------------------------------------------------------------------------------------------------
#                                       PIPELINE BICONNECTED COMPONENTS
#--------------------------------------------------------------------------------------------------------------------------------

"""
    biconnected_components(EV)

Compute sets of 2-cells on the same 3-cells biconnected components.

The method evaluate the 2-cells ``σ_i`` wich lies on the same 3-cells biconnected
component ``χ_j`` and gives back an array that contains an array of ``σ_i`` for each ``j``.
Do note that if the cochain EV contains more copies of the same 2-cell then
it will be considered like a 3-cell.

Intermediate part of Planar Arrangement's algorithmic pipeline.

See also: [`Lar.planar_arrangement`](@ref) for the complete pipeline.

---

# Examples
```jldoctest
julia> EV = SparseArrays.sparse(Array{Int8, 2}([
    [1 1 0 0 0 0] #1 -> 1,2  |
    [1 0 1 0 0 0] #2 -> 1,3  |
    [1 0 0 1 0 0] #3 -> 1,4   |
    [1 0 0 0 1 0] #4 -> 1,5   |
    [1 0 0 0 0 1] #5 -> 1,6
    [0 1 1 0 0 0] #6 -> 2,3  |
    [0 0 0 1 1 0] #7 -> 4,5   |
    ]));
julia> Lar.Arrangement.biconnected_components(EV)
2-element Array{Array{Int64,1},1}:
 [2, 6, 1]
 [4, 7, 3]
```

```jldoctest
julia> EV = SparseArrays.sparse(Array{Int8, 2}([
    [1 1 0] #1 -> 1,2  |
    [1 1 0] #2 -> 1,2  |
    [1 0 1] #3 -> 1,2
    ]));
julia> Lar.Arrangement.biconnected_components(EV)
1-element Array{Array{Int64,1},1}:
 [2, 1]
```
"""
function biconnected_components(EV::Lar.ChainOp)::Array{Array{Int, 1}, 1}
    ps = Array{Tuple{Int, Int, Int}, 1}()
    es = Array{Tuple{Int, Int}, 1}()
    todel = Array{Int, 1}()
    visited = Array{Int, 1}()
    bicon_comps = Array{Array{Int, 1}, 1}()
    hivtx = 1
    
    function an_edge(point)
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



#--------------------------------------------------------------------------------------------------------------------------------
#                                       PIPELINE PART 2
#--------------------------------------------------------------------------------------------------------------------------------

"""
    get_external_cycle(V, EV, FE)

Evaluates the index of the external 3D-cell.

This method looks for and retrieve the external cycle of the model `(V, EV, FE)`.
If the cell does not exist then it return a `Nothing` element.

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest
# ``K^4`` graph
julia> V = [0.0 0.0; 4.0 0.0; 2.0 3.0; 2.0 1.5];
julia> EV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0] #1 -> 1,2
		[0 1 1 0] #2 -> 2,3
		[1 0 1 0] #3 -> 3,1
		[1 0 0 1] #4 -> 1,2
		[0 1 0 1] #5 -> 2,3
		[0 0 1 1] #6 -> 3,1
	]));
julia> FE = SparseArrays.sparse(Array{Int8, 2}([
		[1 0 0 1 1 0] #1 -> 1,4,5
		[0 1 0 0 1 1] #2 -> 2,5,6
		[0 0 1 1 0 1] #3 -> 3,4,6
		[1 1 1 0 0 0] #4 -> 1,2,3 External
	]));

julia> Lar.Arrangement.get_external_cycle(V, EV, FE)
4

julia> typeof(Lar.Arrangement.get_external_cycle(V, EV, FV[1:3,:]))
Nothing
```
"""
function get_external_cycle(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)
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
    pre_containment_test(bboxes)

More goes HERE...

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest

```
"""
function pre_containment_test(bboxes)
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
    prune_containment_graph(n, V, EVs, shells, graph)

More goes HERE...

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest

```
"""
function prune_containment_graph(n, V, EVs, shells, graph)
    
    for i in 1:n
        an_edge = shells[i].nzind[1]
        origin_index = EVs[i][an_edge, :].nzind[1]
        origin = V[origin_index, :]
 
        for j in 1:n
            if i != j
                if graph[i, j] == 1
                    shell_edge_indexes = shells[j].nzind
                    ev = EVs[j][shell_edge_indexes, :]

                    if !Lar.point_in_face(origin, V, ev)
                        graph[i, j] = 0
                    end
                end
             end
         end

     end
     return graph
end


"""
    transitive_reduction!(graph)

More goes HERE...

See also: [`Lar.Arrangement.componentgraph`](@ref).

---

# Examples
```jldoctest

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


"""
    cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)

Cells composing for the Topological Gift Wrapping algorithm.

This is the online part of the TGW algorithm.

See also: [`Lar.Arrangement.planar_arrangement_2`](@ref).

---

# Examples
```jldoctest

```
"""
function cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    function bboxes(V::Lar.Points, indexes::Lar.ChainOp)
        boxes = Array{Tuple{Any, Any}}(undef, indexes.n)
        for i in 1:indexes.n
            v_inds = indexes[:, i].nzind
            boxes[i] = Lar.bbox(V[v_inds, :])
        end
        boxes
    end
    # initialization
    sums = Array{Tuple{Int, Int, Int}}(undef, 0);
	# assembling child components with father components  
    for father in 1:n
        if sum(containment_graph[:, father]) > 0
            father_bboxes = bboxes(V, abs.(EVs[father]')*abs.(boundaries[father]'))
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


"""
    componentgraph(V, copEV, bicon_comps)

Topological Gift Wrapping algorithm on 2D skeletons.

This is the offline part of the TGW algorithm. It takes in input a model and its
biconnected components mapping and evaluates usefull informations:
 - number of biconnected components
 - ?
 - The 1-cells structure (UNMODIFIED). <----------------------------- Could be removed?
 - Association between non-dangling 2-cells and their orientation.
 - Association between 3-cells and 2-cells (with orientation).
 - Association between 3-cells and their orientation.
 - Shell bounding boxes


See also: [`Lar.Arrangement.planar_arrangement_2`](@ref) for the TGW.

---

# Examples
```jldoctest

```
"""
function componentgraph(V, copEV, bicon_comps)

	# arrangement of isolated components
	n = size(bicon_comps, 1)
	shells = Array{Lar.Chain, 1}(undef, n)
	boundaries = Array{Lar.ChainOp, 1}(undef, n)
	EVs = Array{Lar.ChainOp, 1}(undef, n)
	# for each component
	for p=1:n
		ev = copEV[sort(bicon_comps[p]), :]
		# computation of 2-cells 
		fe = Lar.Arrangement.minimal_2cycles(V, ev) 
		# exterior cycle
		global shell_num = Lar.Arrangement.get_external_cycle(V, ev, fe)
		# decompose each fe (co-boundary local to component)
		EVs[p] = ev 
		global tokeep = setdiff(1:fe.m, shell_num)
		boundaries[p] = fe[tokeep, :]
		shells[p] = fe[shell_num, :]
	end
	
#	@show shell_num;
#	@show [SparseArrays.findnz(EVs[k]) for k=1:length(EVs)];
#	@show tokeep;
#	@show [SparseArrays.findnz(boundaries[k]) for k=1:length(boundaries)];
#	@show [SparseArrays.findnz(shells[k]) for k=1:length(shells)];

	# computation of bounding boxes of isolated components
	shell_bboxes = []
	for i in 1:n
		vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
		@show vs_indexes
		push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
		@show shell_bboxes
	end
	# computation and reduction of containment graph
	containment_graph = Lar.Arrangement.pre_containment_test(shell_bboxes)
	@show 1,Matrix(containment_graph)
	containment_graph = Lar.Arrangement.prune_containment_graph(n, V, EVs, shells, containment_graph)
	@show 2,Matrix(containment_graph)
	Lar.Arrangement.transitive_reduction!(containment_graph) 
	@show 3,Matrix(containment_graph)
	return n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
end


"""
    cleandecomposition(V, copEV, sigma)

!!! NOT WORKING

This function clears the model `(V, copEV)` from all edges outside ``σ``.

This function takes a model `(V, copEV)` and a `Lar.Chain` ``σ`` and
gives back the dropped model and and a vector `todel` of the 2-cells to drop.

---

# Examples
```jldoctest
julia> V = [0.0 0.0; 2.0 0.0; 4.0 0.0; 1.0 1.5; 3.0 1.5; 2.0 3.0; 2.0 -3.; -2. 3.0; 6.0 3.0];

julia> EV = SparseArrays.sparse(Array{Int8, 2}([
    [1 0 1 0 0 0 0 0 0] #1 -> 1,3
    [1 0 0 0 0 1 0 0 0] #2 -> 1,6
    [0 0 1 0 0 1 0 0 0] #3 -> 3,6
    [0 1 0 1 0 0 0 0 0] #4 -> 2,4
    [0 1 0 0 1 0 0 0 0] #5 -> 2,5
    [0 0 0 1 1 0 0 0 0] #6 -> 4,5
    [0 0 0 0 0 0 1 1 0] #7 -> 7,8
    [0 0 0 0 0 0 1 0 1] #8 -> 7,9
    [0 0 0 0 0 0 0 1 1] #9 -> 8,9
    ]));

julia> σ = SparseArrays.sparse([0; 0; 0; 1; 1; 1; 0; 0; 0]);

todel, V, EV = Lar.Arrangement.cleandecomposition(V, EV, convert(Lar.Chain, σ))

Plasm.view(convert(Lar.Points, V'), Lar.cop2lar(EV));
```
"""
function cleandecomposition(V::Lar.Points, copEV::Lar.ChainOp, sigma::Lar.Chain)
    # Deletes edges outside sigma area
	todel = []
	new_edges = []
	map(i->new_edges=union(new_edges, edge_map[i]), sigma.nzind) # ???
	ev = copEV[new_edges, :]
	for e in 1:copEV.m
		if !(e in new_edges)  # ???  remove if ???

			vidxs = copEV[e, :].nzind
			v1, v2 = map(i->V[vidxs[i], :], [1,2])
			centroid = .5*(v1 + v2)
			
			if ! Lar.point_in_face(centroid, V, ev) 
				push!(todel, e)
			end
		end
	end
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
	V, copEV = Lar.delete_edges(todel, V, copEV)
    
    # biconnected components
    bicon_comps = Lar.Arrangement.biconnected_components(copEV) # -> arrays of edge indices
    @show bicon_comps,0
    if isempty(bicon_comps)
        println("No biconnected components found.")
        if (return_edge_map)
            return (nothing, nothing, nothing, nothing)
        else
            return (nothing, nothing, nothing)
        end
    end

	#remove dangling edges    
    edges = sort(union(bicon_comps...))
    todel = sort(setdiff(collect(1:size(copEV,1)), edges))
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
    return todel, V, copEV
end


#--------------------------------------------------------------------------------------------------------------------------------
#								MAIN PIPELINE
#--------------------------------------------------------------------------------------------------------------------------------

"""
	function planar_arrangement_1(V, copEV[, sigma[, return_edge_map[, multiproc]]])

First part of arrangement's algorithmic pipeline.

This function computes the pairwise intersection between each edge of a given 2D cellular complex 1-skeleton.
The computation is speeded up via the evaluation of the Spatial Index. See [`Lar.spaceindex`](@ref). 

See also: [`Lar.planar_arrangement`](@ref) for the complete pipeline.

---

# Examples
```jldoctest
julia> using Plasm

julia> EV = [[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [5, 7], [6, 8]];

julia> V = [
           0.0 0.5 0.0 0.5 0.3 1.0 0.3 1.0;
           0.0 0.0 1.0 1.0 0.5 0.5 1.0 1.0
       ];

julia> W = convert(Lar.Points, V'); # Infering type for W = V'

julia> cop_EV = Lar.coboundary_0(EV::Lar.Cells);

julia> W1, copEV1 = Lar.planar_arrangement_1(W::Lar.Points, cop_EV::Lar.ChainOp)
([0.0 0.0; 0.5 0.0; … ; 1.0 0.5; 1.0 1.0], 
  [1 ,  1]  =  1
  [4 ,  1]  =  1
  [1 ,  2]  =  1
  ⋮
  [11,  8]  =  1
  [9 ,  9]  =  1
  [11,  9]  =  1)

julia> EV1 = Lar.cop2lar(copEV1);

julia> V1 = convert(Lar.Points, W1');

julia> Plasm.view(Plasm.numbering(0.1)((V,[[[k] for k=1:size(V,2)], EV])));

julia> Plasm.view(Plasm.numbering(0.1)((V1,[[[k] for k=1:size(V1,2)], EV1])));
```
"""
function planar_arrangement_1(V::Lar.Points, copEV::Lar.ChainOp, 
		sigma::Lar.Chain=spzeros(Int8, 0), 
		return_edge_map::Bool=false, 
		multiproc::Bool=false)

	# data structures initialization
	edgenum = size(copEV, 1)
	edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
	rV = Lar.Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	finalcells_num = 0

	# spaceindex computation
	model = (convert(Lar.Points,V'),Lar.cop2lar(copEV))
	bigPI = Lar.spaceindex(model::Lar.LAR)

	# multiprocessing of edge fragmentation
	if (multiproc == true)
		in_chan = Distributed.RemoteChannel(()->Channel{Int64}(0))
		out_chan = Distributed.RemoteChannel(()->Channel{Tuple}(0))
		ordered_dict = SortedDict{Int64,Tuple}()
		@async begin
			for i in 1:edgenum
				put!(in_chan,i)
			end
			for p in distributed.workers()
				put!(in_chan,-1)
			end
		end
		for p in distributed.workers()
			@async Base.remote_do(frag_edge_channel, p, in_chan, out_chan, V, copEV, bigPI)
		end
		for i in 1:edgenum
			frag_done_job = take!(out_chan)
			ordered_dict[frag_done_job[1]] = frag_done_job[2]
		end
		for (dkey, dval) in ordered_dict
			i = dkey
			v, ev = dval
			newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
			edge_map[i] = newedges_nums
			finalcells_num += size(ev, 1)
			rV, rEV = Lar.skel_merge(rV, rEV, v, ev)
		end
	else 
	# sequential (iterative) processing of edge fragmentation 
		for i in 1:edgenum
			v, ev = Lar.Arrangement.frag_edge(V, copEV, i, bigPI)
			newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
			edge_map[i] = newedges_nums
			finalcells_num += size(ev, 1)
			rV, rEV = Lar.skel_merge(rV, rEV, v, ev) # block diagonal ...
		end
	end
	# merging of close vertices and edges (2D congruence)
	V, copEV = rV, rEV
    if return_edge_map
        V, copEV, edge_map = Lar.Arrangement.merge_vertices!(V, copEV, edge_map)
        return V, copEV, edge_map
    end
	V, copEV = Lar.Arrangement.merge_vertices!(V, copEV, edge_map)
	return V, copEV
end


"""
	function planar_arrangement_2(V, copEV, bicon_comps[, sigma[, return_edge_map[, multiproc]]])

Second part of arrangement's algorithmic pipeline.

This function is the complete Topological Gift Wrapping (TGW) algorithm that is firstly
locally used in order to decompose the 2-cells and then globally to generate the 3-cells
of the arrangement of the ambient space ``E^3``.

During this process each dangling 2-cell is removed.
Do note that the isolated 1-cells are not removed by this procedure.

See also: [`Lar.planar_arrangement`](@ref) for the complete pipeline.

---

# Examples
```jldoctest
julia> V = [0.0 0.0; 2.0 0.0; 4.0 0.0; 1.0 1.5; 3.0 1.5; 2.0 3.0; 3.0 3.0];
julia> W = convert(Lar.Points, V');
julia> EV = SparseArrays.sparse(Array{Int8, 2}([
    [1 1 0 0 0 0 0] #1 -> 1,2
    [0 1 1 0 0 0 0] #2 -> 2,3
    [1 0 0 1 0 0 0] #3 -> 1,4
    [0 0 0 1 0 1 0] #4 -> 4,6
    [0 0 1 0 1 0 0] #5 -> 3,5
    [0 0 0 0 1 1 0] #6 -> 5,6
    [0 1 0 1 0 0 0] #7 -> 2,4
    [0 1 0 0 1 0 0] #8 -> 2,5
    [0 0 0 1 1 0 0] #9 -> 4,5
    [0 0 0 0 0 1 1]
    ]));

julia> Plasm.view(Plasm.numbering(0.5)((W,[[[k] for k=1:size(W,2)], Lar.cop2lar(EV)])));

julia> bicon_comps = Lar.Arrangement.biconnected_components(EV);

julia> V, EV, FE = Lar.Arrangement.planar_arrangement_2(V, EV, bicon_comps);

julia> Plasm.view( Plasm.numbering1(0.5)((V, EV, FE)) );
```
"""
function planar_arrangement_2(V::Lar.Points, copEV::Lar.ChainOp, bicon_comps::Array{Array{Int, 1}, 1}, 
		sigma::Lar.Chain=spzeros(Int8, 0), 
		return_edge_map::Bool=false, 
		multiproc::Bool=false)

	# Topological Gift Wrapping
	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes = 
		Lar.Arrangement.componentgraph(V, copEV, bicon_comps)
	@show containment_graph
	# only in the context of 3D arrangement
	if sigma.n > 0 ## number of columns of sigma > 0
		todel, V, copEV = Lar.Arrangement.cleandecomposition(V, copEV, sigma)
		V, copEV = Lar.delete_edges(todel, V, copEV) # <----------------------non mi è chiaro perchè non è dentro la fun
	end
	# final shell poset aggregation and FE output
	copEV, FE = Lar.Arrangement.cell_merging(
		n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
	if (return_edge_map)
		return V, copEV, FE, edge_map
	end
	return V, copEV, FE
end 
	

"""
    planar_arrangement(V::Points, copEV::ChainOp
        [, sigma::Chain[, return_edge_map::Bool[, multiproc::Bool]]])

Compute the arrangement on the given cellular complex 1-skeleton in 2D.
Whole arrangement's algorithmic pipeline. 

A cellular complex is arranged when the intersection of every possible pair of cell 
of the complex is empty and the union of all the cells is the whole Euclidean space.
The method with no additional arguments specified gives back the full arranged complex
`V`, `EV` and `FE`.

See also:
  - [`Lar.Arrangement.planar_arrangement_1`](@ref)
  - [`Lar.Arrangement.biconnected_components`](@ref)
  - [`Lar.Arrangement.planar_arrangement_2`](@ref)

## Additional arguments:
  - `sigma::Chain`: if specified, `planar_arrangement` will delete from
    the output every edge and face outside this cell. Defaults to an empty cell.
  - `return_edge_map::Bool`: makes the function return also an `edge_map` which
    maps the edges of the imput to the one of the output. Defaults to `false`.
  - `multiproc::Bool`: Runs the computation in parallel mode. Defaults to `false`.

Many examples could be found in Lar Documentation.
"""
function planar_arrangement( V::Lar.Points, copEV::Lar.ChainOp, 
		sigma::Lar.Chain=spzeros(Int8, 0), 
		return_edge_map::Bool=false, 
		multiproc::Bool=false)

	# edge subdivision
	V, copEV = Lar.Arrangement.planar_arrangement_1(V::Lar.Points, copEV::Lar.ChainOp)
	# biconnected components
	bicon_comps = Lar.Arrangement.biconnected_components(copEV::Lar.ChainOp)
	# 2-complex and containment graph
	V, copEV, copFE = Lar.Arrangement.planar_arrangement_2(V::Lar.Points, copEV::Lar.ChainOp, bicon_comps)
	return V, copEV, copFE
end