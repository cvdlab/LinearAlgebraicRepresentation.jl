using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, SparseArrays


#-----------------------------------------------------------------

"""
    planar_arrangement(V::Points, EV::ChainOp, [sigma::Chain], [return_edge_map::Bool], [multiproc::Bool])

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
function planar_arrangement(
        V::Lar.Points, 	
        copEV::Lar.ChainOp, 
        sigma::Lar.Chain=spzeros(Int8, 0), 
        return_edge_map::Bool=false, 
        multiproc::Bool=false)
    
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
            for p in distributed.workers()
                put!(in_chan,-1)
            end
        end
        
        for p in distributed.workers()
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
    
    V, copEV = rV, rEV

    V, copEV = Lar.Arrangement.merge_vertices!(V, copEV, edge_map)
    
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


function cop2lar(cop::Lar.ChainOp)::Lar.Cells
	[findnz(cop[k,:])[1] for k=1:size(cop,1)]
end


#-----------------------------------------------------------------
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, SparseArrays

function randomcuboids(n,scale=1.)
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		sizes = rand(Float64, 2)
		V,(_,EV,_) = Lar.cuboid(corner,true,corner+sizes)
		center = (corner + corner+sizes)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end

V,EV = randomcuboids(10, .5)
V = Plasm.normalize(V,flag=true)
Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])

#model = lar2tria2lar(U,[copEU, copFUE])  ## To finish in utilities.jl

V = convert(Lar.Points, V')
TVs = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewlarcolor(V::Lar.Points, TVs::Array{Lar.Cells})
#-----------------------------------------------------------------
