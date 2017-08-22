include("./utilities.jl")
using NearestNeighbors
include("./minimal_cycles.jl")

function frag_edge(V::Verts, EV::Cells, edgenum::Int)
    alphas = Dict{Float64, Int}()
    edge = EV[edgenum, :]
    for i in 1:size(EV, 1)
        if i != edgenum
            intersection = intersect_edges(V, edge, EV[i, :])
            for (point, alpha) in intersection
                V = [V; point]
                alphas[alpha] = size(V, 1)
            end
        end
    end

    alphas[0.0], alphas[1.0] = edge.nzind

    alphas_keys = sort(collect(keys(alphas)))
    cells_num = length(alphas_keys)-1
    verts_num = size(V, 1)
    EV = spzeros(Int8, cells_num, verts_num)

    for i in 1:cells_num
        EV[i, alphas[alphas_keys[i]]] = 1
        EV[i, alphas[alphas_keys[i+1]]] = 1
    end

    V, EV
end
function intersect_edges(V::Verts, edge1::Cell, edge2::Cell)
    x1, y1, x2, y2 = vcat(map(c->V[c, :], edge1.nzind)...)
    x3, y3, x4, y4 = vcat(map(c->V[c, :], edge2.nzind)...)
    ret = Array{Tuple{Verts, Float64}, 1}()
    denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
    a = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
    b = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
    
    if 0 <= a <= 1 && 0 <= b <= 1
        p = [(x1 + a*(x2-x1))  (y1 + a*(y2-y1))]
        push!(ret, (p, a))
    elseif isnan(a) && isnan(b) 
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
        
    end
    return ret
end
function merge_vertices(V::Verts, EV::Cells, err=1e-4)
    kdtree = KDTree(V')
    tocheck = collect(size(V,1):-1:1)
    todelete = Array{Int64, 1}()
    while !isempty(tocheck)
        vi = pop!(tocheck)
        if !(vi in todelete)
            nearvs = inrange(kdtree, V[vi, :], err)
            for vj in nearvs
                if vj != vi
                    push!(todelete, vj)
                    EV[:,vi] = EV[:, vi] + EV[:, vj]
                end
            end
        end
    end
    
    tokeep = setdiff(collect(1:size(V,1)), todelete)
    EV = EV[:, tokeep]
    V = V[tokeep, :]
    
    tokeep = Array{Int64, 1}()
    cells = [Set(EV[i, :].nzind) for i in size(EV,1):-1:1]
    i = 0
    while !isempty(cells)
        i += 1
        c = pop!(cells)
        if !(length(c) != 2 || c in cells)
            push!(tokeep, i)
        end
    end
    EV = EV[tokeep, :]
    
    V,EV
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
    FV = abs(FE)*EV
    vs = sparsevec(mapslices(sum, abs(EV), 1)).nzind
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
                    contains = false
                    hits = 0
                    visited_verts = []
                    shell_edge_indexes = shells[j].nzind
                    ev = EVs[j][shell_edge_indexes, :]
                    
                    for edge in 1:ev.m
                        a_id, b_id = ev[edge, :].nzind
                        a = V[a_id, :]
                        b = V[b_id, :]
                        v = b - a
                        alpha = (origin[2] - a[2]) / v[2]
                        if 0 <= alpha <= 1
                            x_int = a[1] + v[1]*alpha
                            if x_int > origin[1]
                                if 0 <= alpha <= 1
                                    hits += 1
                                else
                                    p = (alpha == 0) ? a : b
                                    if !(p in visited_verts)
                                        hits += 1
                                        push!(visited_verts, p)
                                    end
                                end
                            end
                        end
                    end
                    
                    contains = hits % 2 == 1
                    
                    if !contains
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
            father_bboxes = bboxes(V, abs(EVs[father]')*abs(boundaries[father]'))
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


function planar_arrangement(V::Verts, EV::Cells)
    EVs = Array{Cells, 1}()
    finalcells_num = 0
    
    edgenum = size(EV, 1)
    for i in 1:edgenum
        V, ev = frag_edge(V, EV, i)
        finalcells_num += size(ev, 1)
        push!(EVs, ev)
        
    end
    EV = spzeros(Int8, finalcells_num, size(V,1))
    newcell_index = 1
    for ev in EVs
        s = size(ev)
        EV[newcell_index:newcell_index+s[1]-1, 1:s[2]] = ev
        newcell_index += s[1]
    end
    
    V, EV = merge_vertices(V, EV)
    
    bicon_comps = biconnected_components(EV)
    
    if isempty(bicon_comps)
        println("No biconnected components found.")
        return
    end
    
    todel = setdiff(collect(1:size(EV,1)), union(bicon_comps...))
    EV = EV[union(bicon_comps...), :]
    
    vertinds = 1:size(EV, 2)
    todel = Array{Int64, 1}()
    for i in vertinds
        if length(EV[:, i].nzind) == 0
            push!(todel, i)
        end
    end
    tokeep = setdiff(vertinds, todel)
    EV = EV[:, tokeep]
    V = V[tokeep, :]
    
    bicon_comps = biconnected_components(EV)
    
    n = size(bicon_comps, 1)
    shells = Array{Cell, 1}(n)
    boundaries = Array{Cells, 1}(n)
    EVs = Array{Cells, 1}(n)
    for p in 1:n
        ev = EV[bicon_comps[p], :]
        fe = minimal_2cycles(V, ev)
        shell_num = get_external_cycle(V, ev, fe)
    
        EVs[p] = ev 
        tokeep = setdiff(1:fe.m, shell_num)
        boundaries[p] = fe[tokeep, :]
        shells[p] = fe[shell_num, :]
    end
    
    shell_bboxes = []
    for i in 1:n
        vs_indexes = (abs(EVs[i]')*abs(shells[i])).nzind
        push!(shell_bboxes, bbox(V[vs_indexes, :]))
    end
    
    containment_graph = pre_containment_test(shell_bboxes)
    containment_graph = prune_containment_graph(n, V, EVs, shells, containment_graph)
    
    transitive_reduction!(containment_graph) 
    
    EV, FE = cell_merging(n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    
    
    

    V, EV, FE
end 
