using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, SparseArrays, DataStructures
using Distributed


#-----------------------------------------------------------------

function splitedges(V,copEV,sigma,return_edge_map,multiproc)
@show multiproc

    edgenum = size(copEV, 1)
    edge_map = Array{Array{Int, 1}, 1}(undef,edgenum)
    rV = Lar.Points(zeros(0, 2))
    rEV = SparseArrays.spzeros(Int8, 0, 0)
    finalcells_num = 0

    if (multiproc == true)
        in_chan = Distributed.RemoteChannel(()->Channel{Int64}(0))
        out_chan = Distributed.RemoteChannel(()->Channel{Tuple}(0))
        ordered_dict = SortedDict{Int64,Tuple}()
        @async begin
            for i in 1:edgenum
                put!(in_chan,i)
            end
            for p in Distributed.workers()
                put!(in_chan,-1)
            end
        end
        for p in Distributed.workers()
            @async Base.remote_do(frag_edge_channel, p, in_chan, out_chan, V, copEV)
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
        for i in 1:edgenum
            v, ev = Lar.Arrangement.frag_edge(V, copEV, i)
            newedges_nums = map(x->x+finalcells_num, collect(1:size(ev, 1)))
            edge_map[i] = newedges_nums
            finalcells_num += size(ev, 1)
            rV, rEV = Lar.skel_merge(rV, rEV, v, ev)
        end
    end
    return V,copEV,sigma,return_edge_map,edge_map,rV, rEV
end

function deleteoutsidesigma(V,copEV,sigma,return_edge_map,edge_map)
    # Deletes edges outside sigma area
    if sigma.n > 0
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
    end
    return V,copEV,sigma,return_edge_map,edge_map
end



function minimal_cycles(angles_fn::Function, verbose=false)

    function _minimal_cycles(V::Lar.Points, 
    ld_bounds::Lar.ChainOp)
    
        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int8, ld_cellsnum)
        dir_marks = zeros(Int8, ld_cellsnum)
        d_bounds = spzeros(Int8, ld_cellsnum, 0)
        
        angles = Array{Array{Int64, 1}, 1}(undef, lld_cellsnum)
        
        function get_seed_cell()
            s = -1
            for i in 1:ld_cellsnum
                if count_marks[i] == 0
                    return i
                elseif count_marks[i] == 1 && s < 0
                    s = i
                end
            end
            return s
        end
        for lld in 1:lld_cellsnum
            as = []
            for ld in ld_bounds[lld, :].nzind
                push!(as, (ld, angles_fn(lld, ld)))
            end
            sort!(as, lt=(a,b)->a[2]<b[2])
            as = map(a->a[1], as)
            angles[lld] = as
        end
        
        function nextprev(lld::Int64, ld::Int64, norp)
            as = angles[lld]
            #ne = findfirst(as, ld)  (findfirst(isequal(v), A), 0)[1]
            ne = (findfirst(isequal(ld), as), 0)[1]  
            while true
                ne += norp
                if ne > length(as)
                    ne = 1
                elseif ne < 1
                    ne = length(as)
                end
        
                if count_marks[as[ne]] < 2
                    break
                end
            end
            as[ne]
        end
        
        
        while (sigma = get_seed_cell()) > 0
        
            if verbose
                print(Int(floor(50 * sum(count_marks) / ld_cellsnum)), "%\r")
            end
        
            c_ld = spzeros(Int8, ld_cellsnum)
            if count_marks[sigma] == 0
                c_ld[sigma] = 1
            else
                c_ld[sigma] = -dir_marks[sigma]
            end
            c_lld = ld_bounds*c_ld
            while c_lld.nzind != []
                corolla = spzeros(Int8, ld_cellsnum)
                for tau in c_lld.nzind
                    b_ld = ld_bounds[tau, :]
                    pivot = intersect(c_ld.nzind, b_ld.nzind)[1]
                    adj = nextprev(tau, pivot, sign(-c_lld[tau]))
                    corolla[adj] = c_ld[pivot]
                    if b_ld[adj] == b_ld[pivot]
                        corolla[adj] *= -1
                    end
                end
                c_ld += corolla
                c_lld = ld_bounds*c_ld
            end
            map(s->count_marks[s] += 1, c_ld.nzind)
            map(s->dir_marks[s] = c_ld[s], c_ld.nzind)
            d_bounds = [d_bounds c_ld]
            
        end
        
        return d_bounds
        
    end

    return _minimal_cycles
end


function minimal_2cycles(V::Lar.Points, EV::Lar.ChainOp)

    function edge_angle(v::Int, e::Int)
        edge = EV[e, :]
        v2 = setdiff(edge.nzind, [v])[1]
        x, y = V[v2, :] - V[v, :]
        return atan(y, x)
    end

    for i in 1:EV.m
        j = min(EV[i,:].nzind...)
        EV[i, j] = -1
    end
    VE = convert(Lar.ChainOp, SparseArrays.transpose(EV))

    EF = minimal_cycles(edge_angle)(V, VE)

    return convert(Lar.ChainOp, SparseArrays.transpose(EF))
end



"""
    test_planar_arrangement(V::Points, EV::ChainOp, [sigma::Chain], [return_edge_map::Bool], [multiproc::Bool])

Compute the arrangement on the given cellular complex 1-skeleton in 2D.

A cellular complex is arranged when the intersection of every possible pair of cell 
of the complex is empty and the union of all the cells is the whole Euclidean space.
The basic method of the function without the `sigma`, `return_edge_map` and `multiproc` arguments 
returns the full arranged complex `V`, `EV` and `FE`.

## Additional arguments:
- `sigma::Chain`: if specified, `planar_arrangement` will delete from the output every edge and face outside this cell. Defaults to an empty cell.
- `return_edge_map::Bool`: makes the function return also an `edge_map` which maps the edges of the imput to the one of the output. Defaults to `false`.
- `multiproc::Bool`: Runs the computation in parallel mode. Defaults to `false`.
"""
function test_planar_arrangement(
	V::Lar.Points, 	
	copEV::Lar.ChainOp; 
	sigma::Lar.Chain=spzeros(Int8, 0), 
	return_edge_map::Bool=false, 
	multiproc::Bool=false)

	V,copEV,sigma,return_edge_map,edge_map,rV, rEV = splitedges(
							V,copEV,sigma,return_edge_map,multiproc)
    V, copEV = rV, rEV
    V, copEV = Lar.Arrangement.merge_vertices!(V, copEV, edge_map)


	V,copEV,sigma,return_edge_map,edge_map = deleteoutsidesigma(
							V,copEV,sigma,return_edge_map,edge_map )
							
    bicon_comps = Lar.Arrangement.biconnected_components(copEV)
    
    if isempty(bicon_comps)
        println("No biconnected components found.")
        if (return_edge_map)
            return (nothing, nothing, nothing, nothing)
        else
            return (nothing, nothing, nothing)
        end
    end
    
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
    
    V, copEV = Lar.delete_edges(todel, V, copEV)
    
    bicon_comps = Lar.Arrangement.biconnected_components(copEV)
    
    n = size(bicon_comps, 1)
    shells = Array{Lar.Chain, 1}(undef, n)
    boundaries = Array{Lar.ChainOp, 1}(undef, n)
    EVs = Array{Lar.ChainOp, 1}(undef, n)
    for p in 1:n
        ev = copEV[sort(bicon_comps[p]), :]
        fe = minimal_2cycles(V, ev)
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
    containment_graph =  Lar.Arrangement.prune_containment_graph(n, V, EVs, shells, containment_graph)
    
     Lar.Arrangement.transitive_reduction!(containment_graph) 
    
    copEV, FE = Lar.Arrangement.cell_merging(
    	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes)
    
    if (return_edge_map)
        return V, copEV, FE, edge_map
    else
        return V, copEV, FE
    end
end 



#-----------------------------------------------------------------

V,EV = Lar.randomcuboids(10, .75)
V = Plasm.normalize(V,flag=true)
Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = test_planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp; 
	multiproc=false)

triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
V = convert(Lar.Points, V')
FVs = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

#-----------------------------------------------------------------
