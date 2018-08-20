function frag_face_channel(in_chan, out_chan, V, EV, FE, sp_idx)
    run_loop = true
    while run_loop
        sigma = take!(in_chan)
        if sigma != -1
            put!(out_chan, frag_face(V, EV, FE, sp_idx, sigma))
        else
            run_loop = false
        end
    end
end
function frag_face(V, EV, FE, sp_idx, sigma)
    vs_num = size(V, 1)

    sigmavs = (abs.(FE[sigma:sigma,:])*abs.(EV))[1,:].nzind 
    sV = V[sigmavs, :]
    sEV = EV[FE[sigma, :].nzind, sigmavs]

    M = submanifold_mapping(sV)
    tV = ([V ones(vs_num)]*M)[:, 1:3]
    
    sV = tV[sigmavs, :]
    
    for i in sp_idx[sigma]
        tmpV, tmpEV = face_int(tV, EV, FE[i, :])
        
        sV, sEV = skel_merge(sV, sEV, tmpV, tmpEV)
    end
    
    sV = sV[:, 1:2]
    

    nV, nEV, nFE = planar_arrangement(sV, sEV, sparsevec(ones(Int8, length(sigmavs))))

    if nV == nothing
        return [], spzeros(Int8, 0,0), spzeros(Int8, 0,0)
    end

    nvsize = size(nV, 1)
    nV = [nV zeros(nvsize) ones(nvsize)]*inv(M)[:, 1:3]
    return nV, nEV, nFE
end

function merge_vertices(V::LinearAlgebraicRepresentation.Points, EV::LinearAlgebraicRepresentation.ChainOp, FE::LinearAlgebraicRepresentation.ChainOp, err=1e-4)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    facenum = size(FE, 1)
    newverts = zeros(Int, vertsnum)
    # KDTree constructor needs an explicit array of Float64
    V = Array{Float64,2}(V)
    kdtree = KDTree(V')

    todelete = []
    
    i = 1
    for vi in 1:vertsnum
        if !(vi in todelete)
            nearvs = LinearAlgebraicRepresentation.inrange(kdtree, V[vi, :], err)
    
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
    
    faces = [[
        map(x->newverts[x], FE[fi, ei] > 0 ? oedges[ei] : reverse(oedges[ei]))
        for ei in FE[fi, :].nzind
    ] for fi in 1:facenum]
    
    
    visited = []
    function filter_fn(face)
    
        verts = []
        map(e->verts = union(verts, collect(e)), face)
        verts = Set(verts)
    
        if !(verts in visited)
            push!(visited, verts)
            return true
        end
        return false
    end
    
    nfaces = filter(filter_fn, faces)
    
    nfacenum = length(nfaces)
    nFE = spzeros(Int8, nfacenum, size(nEV, 1))
    
    for fi in 1:nfacenum
        for edge in nfaces[fi]
            ei = etuple2idx[Tuple{Int, Int}(sort(collect(edge)))]
            nFE[fi, ei] = sign(edge[2] - edge[1])
        end
    end
    
    return LinearAlgebraicRepresentation.Points(nV), nEV, nFE
end


"""
    spatial_arrangement(V::Points, EV::ChainOp, FE::ChainOp; [multiproc::Bool])

Compute the arrangement on the given cellular complex 2-skeleton in 3D.

A cellular complex is arranged when the intersection of every possible pair of cell 
of the complex is empty and the union of all the cells is the whole Euclidean space.
The function returns the full arranged complex as a list of vertices V and a chain of borders EV, FE, CF.

## Additional arguments:
- `multiproc::Bool`: Runs the computation in parallel mode. Defaults to `false`.
"""
function spatial_arrangement(V::LinearAlgebraicRepresentation.Points, EV::LinearAlgebraicRepresentation.ChainOp, FE::LinearAlgebraicRepresentation.ChainOp, multiproc::Bool=false)

    fs_num = size(FE, 1)
    sp_idx = spatial_index(V, EV, FE)

    global rV = LinearAlgebraicRepresentation.Points(undef, 0,3)
    global rEV = SparseArrays.spzeros(Int8,0,0)
    global rFE = SparseArrays.spzeros(Int8,0,0)

    if (multiproc == true)
        in_chan = Distributed.RemoteChannel(()->Channel{Int64}(0))
        out_chan = Distributed.RemoteChannel(()->Channel{Tuple}(0))
        
        @async begin
            for sigma in 1:fs_num
                put!(in_chan, sigma)
            end
            for p in Distributed.workers()
                put!(in_chan, -1)
            end
        end
        
        for p in Distributed.workers()
            @async Base.remote_do(
                frag_face_channel, p, in_chan, out_chan, V, EV, FE, sp_idx)
        end
        
        for sigma in 1:fs_num
            rV, rEV, rFE = skel_merge(rV, rEV, rFE, take!(out_chan)...)
        end
        
    else
        for sigma in 1:fs_num
            # print(sigma, "/", fs_num, "\r")
            nV, nEV, nFE = frag_face(V, EV, FE, sp_idx, sigma)
            rV, rEV, rFE = skel_merge(rV, rEV, rFE, nV, nEV, nFE)
        end
        
    end

    rV, rEV, rFE = merge_vertices(rV, rEV, rFE)
    
    rCF = minimal_3cycles(rV, rEV, rFE)

    return rV, rEV, rFE, rCF
end

