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
model2d = V,EV
Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

#-----------------------------------------------------------------

V,EV
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
Z = convert(Lar.Points, V')
cop_EZ = convert(Lar.ChainOp, cop_EV)
U, copEU, copFUE = planar_arrangement(Z::Lar.Points, cop_EZ::Lar.ChainOp)
W = convert(Lar.Points, U')
EW, FE = cop2lar(copEU), cop2lar(copFUE)
@show W;
@show EW;
@show FE;



open("/tmp/test.txt", "w") do f
    write(f, Lar.lar2obj2D(U,[copEU, copFUE]))
end
model = Lar.obj2lar2D("/tmp/test.txt")
U = convert(Lar.Points, model[1]')

Plasm.view(U,model[2][2])



triangulated_faces = Lar.triangulate2D(model[1], [copEU, copFUE])
TVs = triangulated_faces
hpcs = [ Plasm.lar2hpc(W,TVs[i]) for i=1:length(TVs) ]


#-----------------------------------------------------------------
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, SparseArrays

W = [0.921424 0.783906 1.0 0.862481 0.647673 0.294177 0.545007 0.6026 0.249105 0.373739 0.643302 0.613358 0.74072 0.728415 0.573467 0.561162 0.563934 0.595171 0.734097 0.729797 0.377269 0.372969 0.377079 0.39063 0.561896 0.414482 0.718715 0.513517 0.0943703 0.0739826 0.269193 0.248805 0.262529 0.234926 0.248807 0.0992996 0.149508 0.0 0.744129 0.596022 0.652336 0.50423 0.655607 0.585718 0.688968 0.619079 0.57067 0.287061 0.397911 0.448338 0.164728 0.517317 0.606115 0.554695 0.322323 0.270902; 0.742412 0.710137 0.407609 0.375334 0.540678 0.665496 0.576929 0.413029 0.537848 0.49384 0.528299 0.443496 0.00864823 0.246637 0.0 0.237989 0.184363 0.239748 0.382039 0.530799 0.371723 0.520483 0.378299 0.372109 0.377061 0.372799 0.530479 0.524547 0.420323 0.357549 0.363545 0.30077 0.365709 0.305278 0.762278 0.849756 0.592567 0.680045 0.572753 0.661789 0.420063 0.509099 0.885838 0.811536 0.854458 0.780156 0.289873 0.419416 0.368783 0.0220506 0.151594 0.173066 0.194584 0.406777 0.125814 0.338006]

EW = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 7], [6, 7], [8, 10], [9, 10], [5, 11], [11, 12], [8, 12], [6, 9], [13, 14], [15, 17], [16, 17], [13, 15], [14, 18], [16, 18], [19, 20], [21, 23], [10, 23], [10, 22], [19, 25], [25, 26], [24, 26], [21, 24], [20, 27], [11, 27], [11, 28], [22, 28], [29, 30], [31, 32], [29, 33], [31, 33], [30, 34], [32, 34], [39, 40], [12, 41], [12, 42], [27, 39], [27, 41], [7, 40], [7, 28], [28, 42], [47, 49], [24, 49], [23, 24], [23, 48], [50, 51], [47, 52], [50, 52], [33, 48], [33, 34], [34, 51], [18, 53], [18, 25], [25, 54], [55, 56], [17, 53], [17, 52], [52, 55], [26, 54], [26, 49], [49, 56], [35, 36], [37, 38], [35, 37], [36, 38], [43, 44], [45, 46], [43, 45], [44, 46]]

FE = Array{Int64,1}[[1, 2, 3, 4], [5, 9, 29, 43], [7, 11, 19, 21, 23, 25, 27, 38, 41, 47, 57, 62], [10, 28, 38, 41], [15, 18, 55, 59], [20, 26, 32, 34, 36, 46, 48, 49, 51, 52, 54, 58, 61, 64], [6, 8, 12, 22, 30, 43], [15, 18, 24, 45, 50, 56, 60, 63], [7, 11, 22, 30, 39, 44], [32, 34, 36, 53], [5, 9, 28, 37, 40, 42], [10, 29, 39, 44], [13, 14, 16, 17, 55, 59], [20, 26, 47], [24, 57, 62], [25, 46, 63], [31, 33, 35, 53], [45, 50, 58, 61, 64], [65, 66, 67, 68], [69, 70, 71, 72]]


Plasm.view(W,EW)





