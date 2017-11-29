function frag_edge_channel(in_chan, out_chan, V, EV)
    run_loop = true
    while run_loop
        edgenum = take!(in_chan)
        if edgenum != -1
            put!(out_chan, (edgenum, frag_edge(V, EV, edgenum)))
        else
            run_loop = false
        end
    end
end
function frag_edge(V::Verts, EV::Cells, edge_idx::Int)
    alphas = Dict{Float64, Int}()
    edge = EV[edge_idx, :]
    verts = V[edge.nzind, :]
    for i in 1:size(EV, 1)
        if i != edge_idx
            intersection = intersect_edges(V, edge, EV[i, :])
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
    ev = spzeros(Int8, edge_num, verts_num)

    for i in 1:edge_num
        ev[i, alphas[alphas_keys[i]]] = 1
        ev[i, alphas[alphas_keys[i+1]]] = 1
    end

    verts, ev
end
function intersect_edges(V::Verts, edge1::Cell, edge2::Cell)
    err = 10e-8

    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Verts, Float64}, 1}()

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
function merge_vertices!(V::Verts, EV::Cells, edge_map, err=1e-4)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    newverts = zeros(Int, vertsnum)
    kdtree = KDTree(V')

    todelete = []
    
    i = 1
    for vi in 1:vertsnum
        if !(vi in todelete)
            nearvs = inrange(kdtree, V[vi, :], err)
    
            newverts[nearvs] = i
    
            nearvs = setdiff(nearvs, vi)
            todelete = union(todelete, nearvs)
    
            i = i + 1
        end
    end
    
    nV = V[setdiff(collect(1:vertsnum), todelete), :]
    
    edges = Array{Tuple{Int, Int}, 1}(edgenum)
    oedges = Array{Tuple{Int, Int}, 1}(edgenum)
    
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
        nEV[ei, collect(nedges[ei])] = 1
        etuple2idx[nedges[ei]] = ei
    end
    
    for i in 1:length(edge_map)
        row = edge_map[i]
        row = map(x->edges[x], row)
        row = filter(t->t[1]!=t[2], row)
        row = map(x->etuple2idx[x], row)
        edge_map[i] = row
    end
    

    return nV, nEV
end
function biconnected_components(EV::Cells)
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
        if i == 0
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
function get_external_cycle(V::Verts, EV::Cells, FE::Cells)
    FV = abs.(FE)*EV
    vs = sparsevec(mapslices(sum, abs.(EV), 1)).nzind
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
            if face_area(V, EV, FE[c, :]) < 0
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
            if i != j && bbox_contains(bboxes[j], bboxes[i])
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

                    if !point_in_face(origin, V, ev)
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
    function bboxes(V::Verts, indexes::Cells)
        boxes = Array{Tuple{Any, Any}}(indexes.n)
        for i in 1:indexes.n
            v_inds = indexes[:, i].nzind
            boxes[i] = bbox(V[v_inds, :])
        end
        boxes
    end
    

    sums = Array{Tuple{Int, Int, Int}}(0);

    for father in 1:n
        if sum(containment_graph[:, father]) > 0
            father_bboxes = bboxes(V, abs.(EVs[father]')*abs.(boundaries[father]'))
            for child in 1:n
                if containment_graph[child, father] > 0
                    child_bbox = shell_bboxes[child]
                    for b in 1:length(father_bboxes)
                        if bbox_contains(father_bboxes[b], child_bbox)
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


function planar_arrangement(V::Verts, EV::Cells, sigma::Cell=spzeros(Int8, 0); multiproc=false)
    edgenum = size(EV, 1)
    edge_map = Array{Array{Int, 1}, 1}(edgenum)
    rV = zeros(0, 2)
    rEV = spzeros(Int8, 0, 0)
    finalcells_num = 0
    

    if (multiproc == true)
        in_chan = RemoteChannel(()->Channel{Int64}(0))
        out_chan = RemoteChannel(()->Channel{Tuple}(0))
        
        ordered_dict = SortedDict{Int64,Tuple}()
        
        @schedule begin
            for i in 1:edgenum
                put!(in_chan,i)
            end
            for p in workers()
                put!(in_chan,-1)
            end
        end
        
        for p in workers()
            @async Base.remote_do(frag_edge_channel, p, in_chan, out_chan, V, EV)
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
            
            rV, rEV = skel_merge(rV, rEV, v, ev)
        end
        
    else
        for i in 1:edgenum
            v, ev = frag_edge(V, EV, i)
        
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
        
            edge_map[i] = newedges_nums
        
            finalcells_num += size(ev, 1)
            rV, rEV = skel_merge(rV, rEV, v, ev)
        end
        
    end
    
    V, EV = rV, rEV

    V, EV = merge_vertices!(V, EV, edge_map)
    
    if sigma.n > 0
        todel = []
        
        new_edges = []
        map(i->new_edges=union(new_edges, edge_map[i]), sigma.nzind)
        ev = EV[new_edges, :]
    
        for e in 1:EV.m
            if !(e in new_edges)
    
                vidxs = EV[e, :].nzind
                v1, v2 = map(i->V[vidxs[i], :], [1,2])
                centroid = .5*(v1 + v2)
                
                if !point_in_face(centroid, V, ev) 
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
        
    
        V, EV = delete_edges(todel, V, EV)
    end
    
    bicon_comps = biconnected_components(EV)
    
    if isempty(bicon_comps)
        println("No biconnected components found.")
        return (nothing, nothing, nothing)
    end
    
    
    edges = sort(union(bicon_comps...))
    todel = sort(setdiff(collect(1:size(EV,1)), edges))
    
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
    
    
    V, EV = delete_edges(todel, V, EV)
    
    
    bicon_comps = biconnected_components(EV)
    
    n = size(bicon_comps, 1)
    shells = Array{Cell, 1}(n)
    boundaries = Array{Cells, 1}(n)
    EVs = Array{Cells, 1}(n)
    for p in 1:n
        ev = EV[sort(bicon_comps[p]), :]
        fe = minimal_2cycles(V, ev)
        shell_num = get_external_cycle(V, ev, fe)
    
        EVs[p] = ev 
        tokeep = setdiff(1:fe.m, shell_num)
        boundaries[p] = fe[tokeep, :]
        shells[p] = fe[shell_num, :]
    end
    
    shell_bboxes = []
    for i in 1:n
        vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
        push!(shell_bboxes, bbox(V[vs_indexes, :]))
    end
    
    containment_graph = pre_containment_test(shell_bboxes)
    containment_graph = prune_containment_graph(n, V, EVs, shells, containment_graph)
    
    transitive_reduction!(containment_graph) 
    
    EV, FE = cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    
    
      

    V, EV, FE, edge_map
end 
