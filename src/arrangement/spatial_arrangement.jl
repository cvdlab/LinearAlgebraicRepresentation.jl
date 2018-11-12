Lar = LinearAlgebraicRepresentation
using Plasm

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
        
        sV, sEV = Lar.skel_merge(sV, sEV, tmpV, tmpEV)
    end
    
    sV = sV[:, 1:2]
    

    nV, nEV, nFE = planar_arrangement(sV, sEV, sparsevec(ones(Int, length(sigmavs))))

    if nV == nothing
        return [], spzeros(Int, 0,0), spzeros(Int, 0,0)
    end

    nvsize = size(nV, 1)
    nV = [nV zeros(nvsize) ones(nvsize)]*inv(M)[:, 1:3]
    return nV, nEV, nFE
end

function merge_vertices(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp, err=Lar.ERR)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    facenum = size(FE, 1)
    newverts = zeros(Int, vertsnum)
    
    # KDTree constructor needs an explicit array of Float64 #############
    V = Array{Float64,2}(V)
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
    
    # edge collapse ######################################################
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
    nEV = spzeros(Int, nedgenum, size(nV, 1))
    
    etuple2idx = Dict{Tuple{Int, Int}, Int}()
    
    for ei in 1:nedgenum
        nEV[ei, collect(nedges[ei])] = 1
        etuple2idx[nedges[ei]] = ei
    end
    
    # face collapse ######################################################
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
    nFE = spzeros(Int, nfacenum, size(nEV, 1))
    
    for fi in 1:nfacenum
        for edge in nfaces[fi]
            ei = etuple2idx[Tuple{Int, Int}(sort(collect(edge)))]
            nFE[fi, ei] = sign(edge[2] - edge[1])
        end
    end
    
    return Lar.Points(nV), nEV, nFE
end


"""
    spatial_arrangement(V::Points, EV::ChainOp, FV::Lar.Cells, FE::ChainOp; [multiproc::Bool])

Compute the arrangement on the given cellular complex 2-skeleton in 3D.

A cellular complex is arranged when the intersection of every possible pair of cell 
of the complex is empty and the union of all the cells is the whole Euclidean space.
The function returns the full arranged complex as a list of vertices V and a chain of borders EV, FE, CF.

## Additional arguments:
- `multiproc::Bool`: Runs the computation in parallel mode. Defaults to `false`.
"""
function spatial_arrangement(V::Lar.Points, EV::Lar.ChainOp, FV::Lar.Cells, FE::Lar.ChainOp, multiproc=false)

    fs_num = size(FE, 1)
    sp_idx = Lar.Arrangement.spatial_index(V, EV, FE)

    rV = Lar.Points(0,3)
    rEV = spzeros(Int,0,0)
    rFE = spzeros(Int,0,0)

    if (multiproc == true)
        in_chan = RemoteChannel(()->Channel{Int64}(0))
        out_chan = RemoteChannel(()->Channel{Tuple}(0))
        
        @schedule begin
            for sigma in 1:fs_num
                put!(in_chan, sigma)
            end
            for p in workers()
                put!(in_chan, -1)
            end
        end
        
        for p in workers()
            @async Base.remote_do(
                frag_face_channel, p, in_chan, out_chan, V, EV, FE, sp_idx)
        end
        
        for sigma in 1:fs_num
            rV, rEV, rFE = Lar.skel_merge(rV, rEV, rFE, take!(out_chan)...)
        end
        
    else
		depot_V = Array{Array{Float64,2},1}(fs_num)
		depot_EV = Array{Lar.ChainOp,1}(fs_num)
		depot_FE = Array{Lar.ChainOp,1}(fs_num)
		#V, EV, FE, space_idx = Lar.preprocessing(V,EV,FV)

		for sigma in 1:fs_num
		   print(sigma, "/", fs_num, "\r")
		   #nV, nEV, nFE = Lar.computefragments(V, EV, FV, FE, space_idx, sigma)
		   nV, nEV, nFE = Lar.Arrangement.frag_face(V, EV, FE, sp_idx, sigma)
		   ##@show nV
		   depot_V[sigma] = nV
		   depot_EV[sigma] = Matrix(nEV)
		   depot_FE[sigma] = Matrix(nFE)
		end

		rV = vcat(depot_V...)
		rEV = Lar.blockdiag(depot_EV...)
		rFE = Lar.blockdiag(depot_FE...)
    end

    #rV, rEV, rFE = Lar.merge_vertices(rV, rEV, rFE)
    rV, rEV, rFE = Lar.Arrangement.merge_vertices(rV, rEV, rFE)
		
    rCF = Lar.Arrangement.minimal_3cycles(rV, rEV, rFE)

    return rV, rEV, rFE, rCF
end

function blockdiag(X::SparseMatrixCSC...)
    num = length(X)
    mX = Int[ size(x, 1) for x in X ]
    nX = Int[ size(x, 2) for x in X ]
    m = sum(mX)
    n = sum(nX)

    Tv = promote_type(map(x->eltype(x.nzval), X)...)
    Ti = isempty(X) ? Int : promote_type(map(x->eltype(x.rowval), X)...)

    colptr = Vector{Ti}(n+1)
    nnzX = Int[ nnz(x) for x in X ]
    nnz_res = sum(nnzX)
    rowval = Vector{Ti}(nnz_res)
    nzval = Vector{Tv}(nnz_res)

    nnz_sofar = 0
    nX_sofar = 0
    mX_sofar = 0
    for i = 1 : num
        colptr[(1 : nX[i] + 1) .+ nX_sofar] = X[i].colptr .+ nnz_sofar
        rowval[(1 : nnzX[i]) .+ nnz_sofar] = X[i].rowval .+ mX_sofar
        nzval[(1 : nnzX[i]) .+ nnz_sofar] = X[i].nzval
        nnz_sofar += nnzX[i]
        nX_sofar += nX[i]
        mX_sofar += mX[i]
    end
    colptr[n+1] = nnz_sofar + 1

    SparseMatrixCSC(m, n, colptr, rowval, nzval)
end
