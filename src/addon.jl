function space_arrangement(V::Points, EV::ChainOp, FE::ChainOp, multiproc::Bool=false)

    fs_num = size(FE, 1)
    sp_idx = LinearAlgebraicRepresentation.Arrangement.spatial_index(V, EV, FE)

    global rV = LinearAlgebraicRepresentation.Points(undef, 0,3)
    global rEV = SparseArrays.spzeros(Int,0,0)
    global rFE = SparseArrays.spzeros(Int,0,0)
    
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

#       for sigma in 1:fs_num
#           # print(sigma, "/", fs_num, "\r")
#           nV, nEV, nFE = LinearAlgebraicRepresentation.Arrangement.frag_face(
#           	V, EV, FE, sp_idx, sigma)
#           a,b,c = LinearAlgebraicRepresentation.skel_merge(
#           	rV, rEV, rFE, nV, nEV, nFE)
#           global rV=a; global rEV=b; global rFE=c
#       end

		depot_V = Array{Array{Float64,2},1}(undef,fs_num)
		depot_EV = Array{ChainOp,1}(undef,fs_num)
		depot_FE = Array{ChainOp,1}(undef,fs_num)
		   for sigma in 1:fs_num
			   #print(sigma, "/", fs_num, "\r")
			   nV, nEV, nFE = Arrangement.frag_face( V, EV, FE, sp_idx, sigma)
			   depot_V[sigma] = nV
			   depot_EV[sigma] = nEV
			   depot_FE[sigma] = nFE
		   end
		rV = vcat(depot_V...)
		rEV = SparseArrays.blockdiag(depot_EV...)
		rFE = SparseArrays.blockdiag(depot_FE...)
    
    end

    rV, rEV, rFE = LinearAlgebraicRepresentation.Arrangement.merge_vertices(rV, rEV, rFE)
    
    rCF = Arrangement.minimal_3cycles(rV, rEV, rFE)

    return rV, rEV, rFE, rCF
end



###  2D triangulation
Lar = LinearAlgebraicRepresentation
""" 
	obj2lar2D(path::AbstractString)::LinearAlgebraicRepresentation.LARmodel

Read a *triangulation* from file, given its `path`. Return a `LARmodel` object
"""
function obj2lar2D(path::AbstractString)::LinearAlgebraicRepresentation.LARmodel
    vs = Array{Float64, 2}(undef, 0, 3)
    edges = Array{Array{Int, 1}, 1}()
    faces = Array{Array{Int, 1}, 1}()

    open(path, "r") do fd
		for line in eachline(fd)
			elems = split(line)
			if length(elems) > 0
				if elems[1] == "v"
					x = parse(Float64, elems[2])
					y = parse(Float64, elems[3])
					z = parse(Float64, elems[4])
					vs = [vs; x y z]
				elseif elems[1] == "f"
					# Ignore the vertex tangents and normals
					_,v1,v2,v3 = map(parse, elems)
					append!(edges, map(sort,[[v1,v2],[v2,v3],[v3,v1]]))
					push!(faces, [v1, v2, v3])
				end
				edges = collect(Set(edges))
			end
		end
	end
    return (vs, [edges,faces])::LinearAlgebraicRepresentation.LARmodel
end


""" 
	lar2obj2D(V::LinearAlgebraicRepresentation.Points, 
			cc::LinearAlgebraicRepresentation.ChainComplex)::String

Produce a *triangulation* from a `LARmodel`. Return a `String` object
"""
function lar2obj2D(V::LinearAlgebraicRepresentation.Points, cc::LinearAlgebraicRepresentation.ChainComplex)::String
    assert(length(cc) == 2)
    copEV, copFE = cc
    V = [V zeros(size(V, 1))]

    obj = ""
    for v in 1:size(V, 1)
        obj = string(obj, "v ", 
        	round(V[v, 1], 6), " ", 
        	round(V[v, 2], 6), " ", 
        	round(V[v, 3], 6), "\n")
    end

    triangulated_faces = triangulate2D(V, cc)

	obj = string(obj, "\n")
	for f in 1:copFE.m
		triangles = triangulated_faces[f]
		for tri in triangles
			t = tri
			#t = copCF[c, f] > 0 ? tri : tri[end:-1:1]
			obj = string(obj, "f ", t[1], " ", t[2], " ", t[3], "\n")
		end
	end

    return obj
end


""" 
	triangulate2D(V::LinearAlgebraicRepresentation.Points, 
			cc::LinearAlgebraicRepresentation.ChainComplex)::Array{Any, 1}

Compute a *CDT* for each face of a `ChainComplex`. Return an `Array` of triangles.
"""
function triangulate2D(V::LinearAlgebraicRepresentation.Points, cc::LinearAlgebraicRepresentation.ChainComplex)::Array{Any, 1}
    copEV, copFE = cc
    triangulated_faces = Array{Any, 1}(copFE.m)
	
    for f in 1:copFE.m       
        edges_idxs = copFE[f, :].nzind
        edge_num = length(edges_idxs)
        edges = Array{Int64,1}[] #zeros(Int64, edge_num, 2)

        fv = LinearAlgebraicRepresentation.buildFV(copEV, copFE[f, :])
        vs = V[fv, :]
        
        for i in 1:length(fv)
        	edge = Int64[0,0]
            edge[1] = fv[i]
            edge[2] = i == length(fv) ? fv[1] : fv[i+1]
            push!(edges,edge::Array{Int64,1})
        end
        edges = hcat(edges...)'
        
        triangulated_faces[f] = Triangle.constrained_triangulation(
        vs, fv, edges, fill(true, edge_num))
        tV = V[:, 1:2]
        
        area = LinearAlgebraicRepresentation.face_area(tV, copEV, copFE[f, :])
        if area < 0 
            for i in 1:length(triangulated_faces[f])
                triangulated_faces[f][i] = triangulated_faces[f][i][end:-1:1]
            end
        end
    end

    return triangulated_faces
end


"""
    blockdiag(A...)
Concatenate matrices block-diagonally. Currently only implemented for sparse matrices.
# Examples
```jldoctest
julia> blockdiag(sparse(2I, 3, 3), sparse(4I, 2, 2))
5×5 SparseMatrixCSC{Int64,Int64} with 5 stored entries:
  [1, 1]  =  2
  [2, 2]  =  2
  [3, 3]  =  2
  [4, 4]  =  4
  [5, 5]  =  4
```
"""
