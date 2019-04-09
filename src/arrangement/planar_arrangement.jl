Lar = LinearAlgebraicRepresentation

function frag_edge_channel(in_chan, out_chan, V, EV, bigPI)
    run_loop = true
    while run_loop
        edgenum = take!(in_chan)
        if edgenum != -1
            put!(out_chan, (edgenum, frag_edge(V, EV, edgenum, bigPI)))
        else
            run_loop = false
        end
    end
end

function frag_edge(V::Lar.Points, EV::Lar.ChainOp, edge_idx::Int, bigPI)
    alphas = Dict{Float64, Int}()
    edge = EV[edge_idx, :]
    verts = V[edge.nzind, :]
    for i in bigPI[edge_idx]
        if i != edge_idx
            intersection = Lar.Arrangement.intersect_edges(
            	V, edge, EV[i, :])
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

function intersect_edges(V::Lar.Points, edge1::Lar.Cell, edge2::Lar.Cell)
    err = 10e-8

    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Lar.Points, Float64}, 1}()

    v1 = [x2-x1, y2-y1];
    v2 = [x4-x3, y4-y3];
    v3 = [x3-x1, y3-y1];

    ang1 = dot(normalize(v1), normalize(v2))
    ang2 = dot(normalize(v1), normalize(v3))
    
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

function merge_vertices!(V::Lar.Points, EV::Lar.ChainOp, edge_map, err=1e-4)
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
    
    for i in 1:length(edge_map)
        row = edge_map[i]
        row = map(x->edges[x], row)
        row = filter(t->t[1]!=t[2], row)
        row = map(x->etuple2idx[x], row)
        edge_map[i] = row
    end
    

    return Lar.Points(nV), nEV
end

function biconnected_components(EV::Lar.ChainOp)
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
    
    bicon_comps
end

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

function cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    function bboxes(V::Lar.Points, indexes::Lar.ChainOp)
        boxes = Array{Tuple{Any, Any}}(undef, indexes.n)
        for i in 1:indexes.n
            v_inds = indexes[:, i].nzind
            boxes[i] = Lar.bbox(V[v_inds, :])
        end
        boxes
    end
    

    sums = Array{Tuple{Int, Int, Int}}(undef, 0);

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

    EV = vcat(EVs...)
    edgenum = size(EV, 1)
    facenum = sum(map(x->size(x,1), boundaries))
    FE = spzeros(Int8, facenum, edgenum)
    shells2 = spzeros(Int8, length(shells), edgenum)
    r_offsets = [1]
    c_offset = 1
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
    
    for (f, r, c) in sums
        FE[r_offsets[f]+r-1, :] += shells2[c, :]
    end
    
    return EV, FE
end


function componentgraph(V, copEV)
	# compute containment graph of components
	bicon_comps = Lar.Arrangement.biconnected_components(copEV)

	n = size(bicon_comps, 1)
	shells = Array{Lar.Chain, 1}(undef, n)
	boundaries = Array{Lar.ChainOp, 1}(undef, n)
	EVs = Array{Lar.ChainOp, 1}(undef, n)
	for p in 1:n
		ev = copEV[sort(bicon_comps[p]), :]
		fe = Lar.Arrangement.minimal_2cycles(V, ev)
		shell_num = Lar.Arrangement.get_external_cycle(
			V, ev, fe)

		EVs[p] = ev 
		tokeep = setdiff(1:fe.m, shell_num)
		boundaries[p] = fe[tokeep, :]
		shells[p] = fe[shell_num, :]
	end

	shell_bboxes = []
	for i in 1:n
		vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
		push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
	end
	containment_graph = Lar.Arrangement.pre_containment_test(shell_bboxes)
	containment_graph = Lar.Arrangement.prune_containment_graph(n, V, EVs, shells, containment_graph)
	Lar.Arrangement.transitive_reduction!(containment_graph) 
	return n, containment_graph, V, EVs, boundaries, shells, shell_bboxes
end

function cleandecomposition(V, copEV, sigma)
    # Deletes edges outside sigma area
	todel = []
	new_edges = []
	map(i->new_edges=union(new_edges, edge_map[i]), sigma.nzind)
	ev = copEV[new_edges, :]
	for e in 1:copEV.m
		if !(e in new_edges)

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




"""
    planar_arrangement(V::Points, copEV::ChainOp, [sigma::Chain], [return_edge_map::Bool], [multiproc::Bool])

CCompute the arrangement on the given cellular complex 1-skeleton in 2D.

A cellular complex is arranged when the intersection of every possible pair of cell 
of the complex is empty and the union of all the cells is the whole Euclidean space.
The basic method of the function without the `sigma`, `return_edge_map` and `multiproc` arguments 
returns the full arranged complex `V`, `EV` and `FE`.

## Additional arguments:
- `sigma::Chain`: if specified, `planar_arrangement` will delete from the output every edge and face outside this cell. Defaults to an empty cell.
- `return_edge_map::Bool`: makes the function return also an `edge_map` which maps the edges of the imput to the one of the output. Defaults to `false`.
- `multiproc::Bool`: Runs the computation in parallel mode. Defaults to `false`.
"""
function planar_arrangement(
        V::Lar.Points, 	
        copEV::Lar.ChainOp, 
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
	
	# merging of close vertices
    V, copEV = rV, rEV
    V, copEV = Lar.Arrangement.merge_vertices!(V, copEV, edge_map)
    
	if sigma.n > 0
		todel, V, copEV = cleandecomposition(V, copEV, sigma)
		V, copEV = Lar.delete_edges(todel, V, copEV)
	end
    
	# compute containment graph of 2D components
    n, containment_graph, V, EVs, boundaries, shells, shell_bboxes = componentgraph(V, copEV)

	# Topological Gift Wrapping ?
    copEV, FE = Lar.Arrangement.cell_merging(
    	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    if (return_edge_map)
        return V, copEV, FE, edge_map
    else
        return V, copEV, FE
    end
end 
