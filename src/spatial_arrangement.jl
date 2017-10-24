
include("./utilities.jl")
include("./planar_arrangement.jl")
include("./dimension_travel.jl")
using NearestNeighbors

function merge_vertices(V::Verts, EV::Cells, FE::Cells, err=1e-4)
    vertsnum = size(V, 1)
    edgenum = size(EV, 1)
    facenum = size(FE, 1)
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
        
        edges[ei] = sort([newverts[v1], newverts[v2]])
        oedges[ei] = sort([v1, v2])
    
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
    

    return nV, nEV, nFE
end


function spatial_arrangement(V::Verts, EV::Cells, FE::Cells)
    vs_num = size(V, 1)
    es_num = size(EV, 1)
    fs_num = size(FE, 1)

    sp_idx = spatial_index(V, EV, FE)

    rV = Verts(0,3)
    rEV = spzeros(Int8,0,0)
    rFE = spzeros(Int8,0,0)

    for sigma in 1:fs_num
        println(sigma, "/", fs_num)

        sigmavs = (abs(FE[sigma:sigma,:])*abs(EV))[1,:].nzind 
        sV = V[sigmavs, :]
        sEV = EV[FE[sigma, :].nzind, sigmavs]

        M = submanifold_mapping(sV[1,:], sV[2,:], sV[3,:])
        tV = ([V ones(vs_num)]*M)[:, 1:3]
        
        sV = tV[sigmavs, :]
        
        for i in sp_idx[sigma]
            tmpV, tmpEV = face_int(tV, EV, FE[i, :])
            
            sV, sEV = skel_merge(sV, sEV, tmpV, tmpEV)
        end
        
        sV = sV[:, 1:2]
        

        nV, nEV, nFE = planar_arrangement(sV, sEV, sparsevec(ones(Int8, length(sigmavs))))

        nvsize = size(nV, 1)
        nV = [nV zeros(nvsize) ones(nvsize)]*inv(M)[:, 1:3]

        rV, rEV, rFE = skel_merge(rV, rEV, rFE, nV, nEV, nFE)
    end

    rV, rEV, rFE = merge_vertices(rV, rEV, rFE)

    rCF = minimal_3cycles(rV, rEV, rFE)

    return rV, rEV, rFE, rCF
end

