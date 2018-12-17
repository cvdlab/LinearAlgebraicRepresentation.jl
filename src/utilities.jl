using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

function make_quotients(verts, faces, edges, err=Lar.ERR)
	vertsnum = size(verts, 2)
	kdtree = NearestNeighbors.KDTree(verts)
	todelete = []
	newverts = zeros(Int, vertsnum)
	i = 1
	for vi in 1:vertsnum
		if !(vi in todelete)
			nearvs = NearestNeighbors.inrange(kdtree, verts[:,vi], err)
			newverts[nearvs] = i
			nearvs = setdiff(nearvs, vi)
			todelete = union(todelete, nearvs)
			i = i + 1
		end
	end
	nverts = verts[:,setdiff(collect(1:vertsnum), todelete)]
	vertdict = DataStructures.OrderedDict()
	for k=1:size(nverts,2)
		key = map(Lar.approxVal(Lar.PRECISION), nverts[:,k])
		vertdict[key] = k
	end

	nedges = Array{Array{Int},1}()
	for (v1,v2) in edges
		(w1,w2) = abs.((v1,v2))
		key1 = map(Lar.approxVal(Lar.PRECISION), verts[:,w1])
		key2 = map(Lar.approxVal(Lar.PRECISION), verts[:,w2])
		push!(nedges, [vertdict[key1], vertdict[key2]])
	end
#	nedges = Set(map(sort,nedges))
#	nedges = sort(collect(nedges), by=x->(x[1],x[2]))
#	nedges = filter( x->x[1]!=x[2], nedges )	
	nedges = sort(nedges, by=x->(x[1],x[2]))
	nedges = sort(collect(Set(nedges)), by=x->(x[1],x[2]))
	nedges = filter( x->x[1]!=x[2], nedges )
	nedges = convert(Array{Array{Int64,1},1}, nedges)	

	nfaces = Array{Array{Int},1}()
	for (k,face) in enumerate(faces)
		nface = Int64[]
		for v in face
			key = map(Lar.approxVal(Lar.PRECISION), nverts[:,v])
			push!(nface, vertdict[key])
		end
		if nface ≠ [] 
			nface = sort(collect(Set(nface)))
			push!(nfaces,nface)
		end
	end
	nfaces = collect(Set(nfaces))
	nfaces = sort(nfaces, lt=lexless)
	nfaces = convert(Lar.Cells, nfaces)
	return nverts, nfaces, nedges
end

function interior_to_f(triangle,f,V,FV,EV,FE)
	v1,v2,v3 = triangle
	u = V[:,v2]-V[:,v1]
	v = V[:,v3]-V[:,v1]
	w = cross(u,v)
	T = eye(4); T[1:3,4] = -V[:,v1]
	R = eye(4); R[1:3,1:3] = [u v w]; R = R'
	mapping = R * T
		
	trianglepoints = [[V[:,v1] V[:,v2] V[:,v3]]; ones(3)'] 
	pts = mapping * trianglepoints
	points2D = [pts[r,:] for r = 1:size(pts,1) 
				if !(all(pts[r,:].==0) || all(pts[r,:].==1) )]
	p2D = hcat(points2D...)'
	checkpoint = 0.4995 .* p2D[:,1] + 0.4995 .* p2D[:,2] + 0.001 .* p2D[:,3]
	
	cellpoints = [V[:,FV[f]]; ones(length(FV[f]))' ]
	points = mapping * cellpoints
	verts2D = [points[r,:] for r = 1:size(points,1) 
				if !(all(points[r,:].==0) || all(points[r,:].==1) )]
	P2D = hcat(verts2D...)'

	vdict = Dict(collect(zip(FV[f], 1:length(FV[f]))))
	celledges = [[vdict[v] for v in EV[e]] for e in FE[f]]
	
	out = Lar.pointInPolygonClassification(P2D,celledges)(checkpoint)
	if out=="p_in" 
		return true 
	else 
		return false 
	end
end


function ordering(triangles,V)
	normals = []
	v1,v2,v3 = triangles[1]
	if v1>v2 v1,v2 = v2,v1 end
	e3 = normalize(V[:,v2]-V[:,v1])
	e1 = normalize(V[:,v3]-V[:,v1])
	e2 = normalize(cross(e1,e3))
	basis = [e1 e2 e3]
	transform = inv(basis)
	
	angles = []
	for (v1,v2,v3) in triangles
		w1 = normalize(V[:,v3]-V[:,v1])
		w2 = transform * w1
		w3 = cross([0,0,1],w2)
		push!(normals,w3)
	end
	for k=1:length(normals)
		angle = atan2(normals[k][2],normals[k][1])
		push!(angles,angle)
	end
	#pairs = sort(collect(zip(angles,1:length(triangles))))
	pairs = sort([(angle,k) for (k,angle) in enumerate(angles)])
	order = [k for (angle,k) in pairs]
	return order
end

function ord(hinge::Int64, bd1::SparseVector{Int64,Int64}, V::Array{Float64,2},
FV::Array{Array{Int64,1},1}, EV::Array{Array{Int64,1},1}, FE::Array{Array{Int64,1},1})
	cells = SparseArrays.findnz(bd1)[1]
	triangles = []

	function area(v1,v2,v3)
		u = V[:,v2]-V[:,v1]
		v = V[:,v3]-V[:,v1]
		out = norm(cross(u,v)) # actually, to be divided by two
		return out
	end
	
	for f in cells
		v1,v2 = EV[hinge]		
		index = findfirst(v -> (area(v1,v2,v)≠0), FV[f])
		v3 = FV[f][index]
		
		# test if [v1,v2,v3] interior to f
		while true
			if interior_to_f([v1, v2, v3],f,V,FV,EV,FE)
				push!(triangles, [v1,v2,v3])
				break
			else
				index = findnext(v -> (area(v1,v2,v)≠0), FV[f], index+1)
				v3 = FV[f][index]
			end
		end
	end
	order = ordering(triangles,V)
	return [cells[index] for index in order]
end


function mynext(cycle, pivot)
	len = length(cycle)
	ind = find(x -> x==pivot, cycle)[1]
	nextIndex = ind==len ? 1 : ind+1
	return cycle[nextIndex][1]
end

function myprev(cycle, pivot)
	len = length(cycle)
	ind = find(x->x==pivot, cycle)[1]
	nextIndex = ind==1 ? len : ind-1
	return cycle[nextIndex][1]
end


function build_copFC(V,FV,EV,copFE)
@show V,FV,EV,SparseArrays.findnz(copFE)
	copEF = copFE'
	FE = [SparseArrays.findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	# Initializations
	m,n = size(copEF)
	marks = zeros(Int,n);
	I = Int64[]; J = Int64[]; W = Int8[]; 
	jcol = 0
	choose(marks) = findfirst(x -> x<2, marks)

	# Main loop (adding one copFC's column stepwise)
	while sum(marks) < 2n
		# select a (d−1)-cell, "seed" of the column extraction
		σ = choose(marks)
		@show σ
		@show marks
		if marks[σ] == 0 
			cd1 = sparsevec([σ], Int8[1], n)
		elseif marks[σ] == 1
			cd1 = sparsevec([σ], Int8[-1], n)
		end
		# compute boundary cd2 of seed cell
		cd2 = copEF * cd1
		# loop until (boundary) cd2 becomes empty
		while nnz(cd2)≠0
			corolla = sparsevec([], Int8[], n)
			# for each “hinge” τ cell
			for τ ∈ (.*)(SparseArrays.findnz(cd2)...)
				#compute the  coboundary
				tau = sparsevec([abs(τ)], [sign(τ)], m)
				bd1 = transpose(transpose(tau) * copEF)
				cells2D = SparseArrays.findnz(bd1)[1]
				# compute the  support
				pivot = intersect(cells2D, SparseArrays.findnz(cd1)[1])[1] 
				# compute the new adj cell
				fan = ord(abs(τ),bd1,V,FV,EV,FE) # ord(pivot,bd1)
				if τ > 0 
					adj = mynext(fan,pivot)
				elseif τ < 0 
					adj = myprev(fan,pivot)
				end
				# orient adj
				if copEF[abs(τ),adj] ≠ copEF[abs(τ),pivot] 
					corolla[adj] = cd1[pivot]
				else
					corolla[adj] = -(cd1[pivot])
				end
			end
			# insert corolla cells in current cd1
			for (k,val) in zip(SparseArrays.findnz(corolla)...) 
				cd1[k] = val 
			end 
			# compute again the boundary of cd1
			cd2 = copEF * cd1
		end
		for σ ∈ SparseArrays.findnz(cd1)[1]
			# update the counters of used cells
			marks[σ] += 1
		end
		# append a new column to [∂d+]
		# copFC += cd1
		rows, vals = SparseArrays.findnz(cd1)
		jcol += 1
		append!(I,rows)
		append!(J,[ jcol for k=1:nnz(cd1) ])
		append!(W,vals)
	end
	copFC = sparse(I,J,W)
	return copFC
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

"""
    bbox(vertices::Points)

The axis aligned bounding box of the provided set of n-dim `vertices`.

The box is returned as the couple of `Points` of the two opposite corners of the box.
"""
function bbox(vertices::Points)
    minimum = mapslices(x->min(x...), vertices, 1)
    maximum = mapslices(x->max(x...), vertices, 1)
    minimum, maximum
end

"""
    bbox_contains(container, contained)

Check if the axis aligned bounding box `container` contains `contained`.

Each input box must be passed as the couple of `Points` standing on the opposite corners of the box.
"""
function bbox_contains(container, contained)
    b1_min, b1_max = container
    b2_min, b2_max = contained
    all(map((i,j,k,l)->i<=j<=k<=l, b1_min, b2_min, b2_max, b1_max))
end

"""
    face_area(V::Points, EV::Cells, face::Cell)

The area of `face` given a geometry `V` and an edge topology `EV`.
"""
function face_area(V::Points, EV::Cells, face::Cell)
    return face_area(V, build_copEV(EV), face)
end

function face_area(V::Points, EV::ChainOp, face::Cell)
    function triangle_area(triangle_points::Points)
        ret = ones(3,3)
        ret[:, 1:2] = triangle_points
        return .5*det(ret)
    end

    area = 0

    fv = buildFV(EV, face)

    verts_num = length(fv)
    v1 = fv[1]

    for i in 2:(verts_num-1)

        v2 = fv[i]
        v3 = fv[i+1]

        area += triangle_area(V[[v1, v2, v3], :])
    end

    return area
end

"""
    skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)

Merge two **1-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)
    V = [V1; V2]
    EV = spzeros(Int, EV1.m + EV2.m, EV1.n + EV2.n)
    EV[1:EV1.m, 1:EV1.n] = EV1
    EV[EV1.m+1:end, EV1.n+1:end] = EV2
    V, EV
end

"""
    skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)

Merge two **2-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)
    FE = spzeros(Int, FE1.m + FE2.m, FE1.n + FE2.n)
    FE[1:FE1.m, 1:FE1.n] = FE1
    FE[FE1.m+1:end, FE1.n+1:end] = FE2
    V, EV = skel_merge(V1, EV1, V2, EV2)
    V, EV, FE
end

"""
    delete_edges(todel, V::Points, EV::ChainOp)

Delete edges and remove unused vertices from a **2-skeleton**.

Loop over the `todel` edge index list and remove the marked edges from `EV`.
The vertices in `V` which remained unconnected after the edge deletion are deleted too.
"""
function delete_edges(todel, V::Points, EV::ChainOp)
    tokeep = setdiff(collect(1:EV.m), todel)
    EV = EV[tokeep, :]
    
    vertinds = 1:EV.n
    todel = Array{Int64, 1}()
    for i in vertinds
        if length(EV[:, i].nzind) == 0
            push!(todel, i)
        end
    end

    tokeep = setdiff(vertinds, todel)
    EV = EV[:, tokeep]
    V = V[tokeep, :]

    return V, EV
end

"""
    buildFV(EV::Cells, face::Cell)

The list of vertex indices that expresses the given `face`.

The returned list is made of the vertex indices ordered following the traversal order to keep a coherent face orientation. 
The edges are need to understand the topology of the face.

In this method the input face must be expressed as a `Cell`(=`SparseVector{Int, Int}`) and the edges as `Cells`.
"""
function buildFV(EV::Cells, face::Cell)
    return buildFV(build_copEV(EV), face)
end

"""
    buildFV(copEV::ChainOp, face::Cell)

The list of vertex indices that expresses the given `face`.

The returned list is made of the vertex indices ordered following the traversal order to keep a coherent face orientation. 
The edges are need to understand the topology of the face.

In this method the input face must be expressed as a `Cell`(=`SparseVector{Int, Int}`) and the edges as `ChainOp`.
"""
function buildFV(copEV::ChainOp, face::Cell)
    startv = -1
    nextv = 0
    edge = 0

    vs = Array{Int64, 1}()

    while startv != nextv
        if startv < 0
            edge = face.nzind[1]
            startv = copEV[edge,:].nzind[face[edge] < 0 ? 2 : 1]
            push!(vs, startv)
        else
            edge = setdiff(intersect(face.nzind, copEV[:, nextv].nzind), edge)[1]
        end
        nextv = copEV[edge,:].nzind[face[edge] < 0 ? 1 : 2]
        push!(vs, nextv)

    end

    return vs[1:end-1]
end

"""
    buildFV(copEV::ChainOp, face::Array{Int, 1})

The list of vertex indices that expresses the given `face`.

The returned list is made of the vertex indices ordered following the traversal order to keep a coherent face orientation. 
The edges are need to understand the topology of the face.

In this method the input face must be expressed as a list of vertex indices and the edges as `ChainOp`.
"""
function buildFV(copEV::ChainOp, face::Array{Int, 1})
    startv = face[1]
    nextv = startv

    vs = []
    visited_edges = []

    while true
        curv = nextv
        push!(vs, curv)

        edge = 0
        for edge in copEV[:, curv].nzind
            nextv = setdiff(copEV[edge, :].nzind, curv)[1]
            if nextv in face && (nextv == startv || !(nextv in vs)) && !(edge in visited_edges)
                break
            end
        end

        push!(visited_edges, edge)

        if nextv == startv
            break
        end
    end

    return vs
end


"""
   build_K(FV::Lar.Cells)::ChainOp

The *characteristic matrix* of type `ChainOp` from 1-cells (edges) to 0-cells (vertices)
"""
function build_K(FV::Lar.Cells)
	I = Int64[]; J = Int64[]; V = Int64[]
	for (i,face) in enumerate(FV)
		for v in face
			push!(I,i)
			push!(J,v)
			push!(V,1)
		end
	end
	kEV = sparse(I,J,V)
end

function build_TE(TV, edge_dict)
	global TE = []
	for tria in TV
		v1,v2,v3 = tria
		e1 = get(edge_dict, [v2,v3], 0)
		if e1 == 0 e1 = get(edge_dict, [v3,v2], 0) end
		e2 = get(edge_dict, [v3,v1], 0)
		if e2 == 0 e2 = get(edge_dict, [v1,v3], 0) end
		e3 = get(edge_dict, [v1,v2], 0)
		if e3 == 0 e3 = get(edge_dict, [v2,v1], 0) end
		te = [e1,e2,e3]
		append!(TE, [te])
	end
	return TE
end

function build_ET(TE)
	edges = [[edge for edge in tria if edge≠0] for tria in TE]
	edges = union(cat(edges))
	global ET = Dict(zip(edges,[Int[] for k=1:length(edges)])) 
	for (k,t) in enumerate(TE)
		for h=1:3
			if t[h]≠0 push!(ET[t[h]], k) end
		end
	end
	return ET
end


function remove_interior_triangles(points2D, edges_list, points_map, 
TV, edge_dict, edges,EV)
	# get interior boundary 1-cycles (s.t. |1-cell coboundary| = 2)
	TE = build_TE(TV, edge_dict)
	ET = build_ET(TE)
	inv_map = Dict(zip(points_map, 1:length(points_map)))
	global triangles = []
	for e ∈ edges
		v1,v2 = EV[e]
		e_coboundary = ET[e]
		if length(e_coboundary) == 1
			t = e_coboundary[1]
			push!(triangles, TV[t])
		elseif length(e_coboundary) == 2
			t1,t2 = e_coboundary
			v = setdiff(Set(TV[t1]),[v1,v2])[1]
			w1,w2,w3 = map(x->inv_map[x],[v1,v2,v])
			point = 0.49*points2D'[:,w1]+0.49*points2D'[:,w2]+0.02*points2D'[:,w3]
			ev = [edges_list[k,:] for k=1:size(edges_list,1)]
			ew = [[inv_map[v1],inv_map[v2]] for (v1,v2) in ev]
			
			if Lar.pointInPolygonClassification(points2D', ew)(point)=="p_in"
				push!(triangles, TV[t1])
			else
				push!(triangles, TV[t2])
			end
		else
			error("More then two coboundary 2-cells of a face edge")
		end
	end
	return triangles
end

function fix_copFE(copFE,kFV,kEV,V,FV,EV,copEF)
println("eccomi 4")
	for f=1:size(copFE,1)
		# compare the numbers of incident edges and vertices (on the f face)
		face_edges = sparsevec( SparseArrays.findnz(copFE[f,:])... )
		face_verts = sparsevec( SparseArrays.findnz(kFV[f,:])... )
		if nnz(face_edges) == nnz(face_verts)
			# standard case: go ahead
			copFE,V,FV,EV = build_kFE(V,FV,EV,false)
			return copFE,V,FV,EV
		elseif nnz(face_edges) > nnz(face_verts)
			# compute "external" edges (may exist iff non-convex f)
			triples = []; extern = []
			for e in SparseArrays.findnz(face_edges)[1]
				vert_edges = sparsevec( SparseArrays.findnz(kEV[e,:])... )
				v1,v2 = SparseArrays.findnz(vert_edges)[1]
				push!(triples, [e,v1,1])
				push!(triples, [e,v2,1])
			end
			I,J,Val = collect(zip(triples...))
			I = [i for i in I]
			J = [j for j in J]
			Val = Int8[v for v in Val]
			m,n = size(kEV)
			copev = sparse(I,J,Val, m,n)
			oddverts = Set([v for v in J if nnz(copev[:,v])>2])
			for e in SparseArrays.findnz(face_edges)[1] 
				v1,v2 = SparseArrays.findnz(kEV[e,:])[1]
				@show e,v1,v2
				if intersect( [v1,v2],oddverts )==[v1,v2]
					@show f,e,v1,v2
					push!(extern, [f,e,v1,v2])
				end
			end

			# insert a new vertex in all external edges
			nverts = size(V,2)
			nedges = length(EV)
			for (k,quadruple) in enumerate(extern)
				f,e,v1,v2 = quadruple
				# insert new vert in V
				last_v = nverts + k
				newvert = 0.5 * V[:,v1] + 0.5 * V[:,v2]
				V = [V newvert]
				# insert new edge
				EV[e] = [v1, last_v]
				push!(EV, [v2, last_v])
				# insert new vertex in incident faces
				incfaces = setdiff(findnz(copEF[e,:])[1], [f])
				for k=1:length(incfaces)
					push!(FV[incfaces[k]],last_v)
				end
			end
		else
			# non-manifold face (wrong case: unfits cell definition)
			@assert numedges<length(verts)  "STOP: non manifold $f face"  
 		end
	end # main loop on faces
	# recompute copFE via kFV and kEV product (no external edges now possible)
println("eccomi 5")
	copFE,V,FV,EV = build_kFE(V,FV,EV,false)
println("eccomi 6")
	return copFE,V,FV,EV
end



function build_kFE(V::Lar.Points, FV::Lar.Cells, EV::Lar.Cells, fix=true::Bool)
	#compute some usigned characteristic matrices
println("eccomi 2")
	kEV = Lar.build_K(EV)
	kFV = Lar.build_K(FV)
	kFE = kFV * kEV'	
	# compute the unsigned operator cobuondary_1 (copFE)
	I = Int64[]; J = Int64[]; W = Int64[]
	for (i,j,v) in zip(SparseArrays.findnz(kFE)...)
		if v == 2
			push!(I,i); push!(J,j); push!(W,1)
		end
	end
	copFE = sparse(I,J,W)
	copEF = copFE'
	if fix==true
println("eccomi 3")
		copFE,V,FV,EV = fix_copFE(copFE,kFV,kEV,V,FV,EV,copEF)
println("eccomi 7")
	end
	return copFE,V,FV,EV
end


"""
   build_copFE(V::Lar.Points, FV::Lar.Cells, EV::Lar.Cells, signed=true)::ChainOp

The signed `ChainOp` from 1-cells (edges) to 2-cells (faces)
"""
function build_copFE(V::Lar.Points, FV::Lar.Cells, EV::Lar.Cells, signed=true)
println("eccomi 1")
	copFE,V,FV,EV = build_kFE(V,FV,EV)
println("eccomi 8")
	if !signed return copFE
	# compute the signed operator cobuondary_1 (signed copFE)
	else
		FE = [SparseArrays.findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
		
		# make dictionary (vi,vj)=>k  (i.e. fk)
		edge_dict = Dict([(EV[k],k) for k=1:size(EV,1)])
		global signedFE = []
		
		# main loop on faces
		for f = 1:size(FE,1)
			fe = []
			# prepare the face triangulation
			points_map = FV[f]
		
			# face mapping
			@show f
			T = Lar.face_mapping(V, FV, f)
			points = [hcat([V[:,v] for v in points_map]...); ones(length(points_map))']
			pts = T * points
			points2D = [pts[r,:] for r = 1:size(pts,1) 
							if !(all(pts[r,:].==0) || all(pts[r,:].==1) )]
			points2D = hcat(points2D...)
		
			edges_list = hcat([EV[e] for e in FE[f]]...)'
			edges = FE[f]
			# get the face triangulation
			TV = Triangle.constrained_triangulation(points2D, points_map, edges_list)
		
			# removal of interior triangles
			tv = remove_interior_triangles(points2D, edges_list, points_map, 
					TV, edge_dict, edges,EV)
		
			# interpretation of triangles to signed edge triples
			function get_signed_edge(v1,v2)
				edge = get(edge_dict, [v1,v2], 0)
				if edge==0 
					edge = get(edge_dict, [v2,v1], 0)*-1
				end
				return edge
			end
		
			for tria in tv
				v1,v2,v3 = tria
				e1 = get_signed_edge(v2,v3)
				e2 = get_signed_edge(v3,v1)
				e3 = get_signed_edge(v1,v2)
				append!(fe,[e1,e2,e3])
			end
			push!(signedFE, union(fe))  #BAH ... !!
		end
		signedFE = [[e for e in face if e≠0] for face in signedFE]
		m,n = size(copFE)
		for i=1:m
			for k=1:length(signedFE[i])
				j = abs(signedFE[i][k])
				val = sign(signedFE[i][k])
				copFE[i,j] = val
			end
		end
		
	end
	return copFE,V,FV,EV
end

"""
   build_copEV(EV::Lar.Cells, signed=true)::ChainOp

The signed (or not) `ChainOp` from 0-cells (vertices) to 1-cells (edges)
"""
function build_copEV(EV::Lar.Cells, signed=true)
	copEV = build_K(EV)
	for e=1:SparseArrays.size(copEV,1)
		h,_ = SparseArrays.findnz(copEV[e,:])[1]
		if signed 
			copEV[e,h] = -1
		else
			copEV[e,h] = 1
		end
	end
	return copEV
end


"""
    vin(vertex, vertices_set)

Checks if `vertex` is one of the vertices inside `vertices_set`
"""
function vin(vertex, vertices_set)
    for v in vertices_set
        if vequals(vertex, v)
            return true
        end
    end
    return false
end

"""
    vequals(v1, v2)

Check the equality between vertex `v1` and vertex `v2`
"""
function vequals(v1, v2)
    err = LinearAlgebraicRepresentation.ERR
    return length(v1) == length(v2) && all(map((x1, x2)->-err < x1-x2 < err, v1, v2))
end

"""
    triangulate(model::LARmodel)

Full constrained Delaunnay triangulation of the given 3-dimensional `LARmodel`
"""
function triangulate(model::LARmodel)
    V, topology = model
    cc = build_cops(topology...)
    return triangulate(V, cc)
end

"""
    triangulate(V::Points, cc::ChainComplex)

Full constrained Delaunnay triangulation of the given 3-dimensional model (given with topology as a `ChainComplex`)
"""
function triangulate(V::Points, cc::ChainComplex)
    copEV, copFE = cc

    triangulated_faces = Array{Any, 1}(copFE.m)

    for f in 1:copFE.m
        if f % 10 == 0
            print(".")
        end
        
        edges_idxs = copFE[f, :].nzind
        edge_num = length(edges_idxs)
        edges = zeros(Int64, edge_num, 2)

        
        fv = buildFV(copEV, copFE[f, :])

        vs = V[fv, :]

        v1 = normalize(vs[2, :] - vs[1, :])
        v2 = [0 0 0]
        v3 = [0 0 0]
        err = LinearAlgebraicRepresentation.ERR
        i = 3
        while -err < norm(v3) < err
            v2 = normalize(vs[i, :] - vs[1, :])
            v3 = cross(v1, v2)
            i = i + 1
        end
        M = reshape([v1; v2; v3], 3, 3)

        vs = (vs*M)[:, 1:2]
        
        for i in 1:length(fv)
            edges[i, 1] = fv[i]
            edges[i, 2] = i == length(fv) ? fv[1] : fv[i+1]
        end
        
        triangulated_faces[f] = Triangle.constrained_triangulation(vs, fv, edges, fill(true, edge_num))

        tV = (V*M)[:, 1:2]
        
        area = face_area(tV, copEV, copFE[f, :])
        if area < 0 
            for i in 1:length(triangulated_faces[f])
                triangulated_faces[f][i] = triangulated_faces[f][i][end:-1:1]
            end
        end
    end

    return triangulated_faces
end

"""
    point_in_face(point, V::Points, copEV::ChainOp)

Check if `point` is inside the area of the face bounded by the edges in `copEV`
"""
function point_in_face(point, V::Points, copEV::ChainOp)

    function pointInPolygonClassification(V,EV)

        function crossingTest(new, old, status, count)
        if status == 0
            status = new
            return status, (count + 0.5)
        else
            if status == old
                return 0, (count + 0.5)
            else
                return 0, (count - 0.5)
            end
        end
        end

        function setTile(box)
        tiles = [[9,1,5],[8,0,4],[10,2,6]]
        b1,b2,b3,b4 = box
        function tileCode(point)
            x,y = point
            code = 0
            if y>b1 code=code|1 end
            if y<b2 code=code|2 end
            if x>b3 code=code|4 end
            if x<b4 code=code|8 end
            return code
        end
        return tileCode
        end

        function pointInPolygonClassification0(pnt)
            x,y = pnt
            xmin,xmax,ymin,ymax = x,x,y,y
            tilecode = setTile([ymax,ymin,xmax,xmin])
            count,status = 0,0

            for k in 1:EV.m
                edge = EV[k,:]
                p1, p2 = V[edge.nzind[1], :], V[edge.nzind[2], :]
                (x1,y1),(x2,y2) = p1,p2
                c1,c2 = tilecode(p1),tilecode(p2)
                c_edge, c_un, c_int = xor(c1, c2), c1|c2, c1&c2
                
                if (c_edge == 0) & (c_un == 0) return "p_on" 
                elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
                elseif c_edge == 3
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 count += 1 end
                elseif c_edge == 15
                    x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
                    if x_int > x count += 1
                    elseif x_int == x return "p_on" end
                elseif (c_edge == 13) & ((c1==4) | (c2==4))
                        status, count = crossingTest(1,2,status,count)
                elseif (c_edge == 14) & ((c1==4) | (c2==4))
                        status, count = crossingTest(2,1,status,count)
                elseif c_edge == 7 count += 1
                elseif c_edge == 11 count = count
                elseif c_edge == 1
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 
                        status, count = crossingTest(1,2,status,count) 
                    end
                elseif c_edge == 2
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 
                        status, count = crossingTest(2,1,status,count) 
                    end
                elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
                elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
                elseif c_edge == 5
                    if (c1==0) | (c2==0) return "p_on"
                    else 
                        status, count = crossingTest(1,2,status,count) 
                    end
                elseif c_edge == 6
                    if (c1==0) | (c2==0) return "p_on"
                    else 
                        status, count = crossingTest(2,1,status,count) 
                    end
                elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
                elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
                end
            end
            
            if (round(count)%2)==1 
                return "p_in"
            else 
                return "p_out"
            end
        end
        return pointInPolygonClassification0
    end
    
    return pointInPolygonClassification(V, copEV)(point) == "p_in"
end

"""
    lar2obj(V::Points, cc::ChainComplex)

Triangulated OBJ string representation of the model passed as input.

Use this function to export LAR models into OBJ

# Example

```julia
	julia> cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
	[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
	[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
	
	julia> cube_2 = Lar.Struct([Lar.t(0,0,0.5), Lar.r(0,0,pi/3), cube_1])
	
	julia> V, FV, EV = Lar.struct2lar(Lar.Struct([ cube_1, cube_2 ]))
	
	julia> V, bases, coboundaries = Lar.chaincomplex(V,FV,EV)
	
	julia> (EV, FV, CV), (copEV, copFE, copCF) = bases, coboundaries

	julia> FV # bases[2]
	18-element Array{Array{Int64,1},1}:
	 [1, 3, 4, 6]            
	 [2, 3, 5, 6]            
	 [7, 8, 9, 10]           
	 [1, 2, 3, 7, 8]         
	 [4, 6, 9, 10, 11, 12]   
	 [5, 6, 11, 12]          
	 [1, 4, 7, 9]            
	 [2, 5, 11, 13]          
	 [2, 8, 10, 11, 13]      
	 [2, 3, 14, 15, 16]      
	 [11, 12, 13, 17]        
	 [11, 12, 13, 18, 19, 20]
	 [2, 3, 13, 17]          
	 [2, 13, 14, 18]         
	 [15, 16, 19, 20]        
	 [3, 6, 12, 15, 19]      
	 [3, 6, 12, 17]          
	 [14, 16, 18, 20]        

	julia> CV # bases[3]
	3-element Array{Array{Int64,1},1}:
	 [2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20]
	 [2, 3, 5, 6, 11, 12, 13, 17]                    
	 [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]    
	 
	julia> copEV # coboundaries[1]
	34×20 SparseMatrixCSC{Int,Int64} with 68 stored entries: ...

	julia> copFE # coboundaries[2]
	18×34 SparseMatrixCSC{Int,Int64} with 80 stored entries: ...
	
	julia> copCF # coboundaries[3]
	4×18 SparseMatrixCSC{Int,Int64} with 36 stored entries: ...
	
	objs = Lar.lar2obj(V'::Lar.Points, [coboundaries...])
			
	open("./two_cubes.obj", "w") do f
    	write(f, objs)
	end


```
"""
function lar2obj(V::Points, cc::ChainComplex)
    copEV, copFE, copCF = cc

    obj = ""
    for v in 1:size(V, 1)
        obj = string(obj, "v ", round(V[v, 1], 6), " ", round(V[v, 2], 6), " ", 
        	round(V[v, 3], 6), "\n")
    end

    print("Triangulating")
    triangulated_faces = triangulate(V, cc[1:2])
    println("DONE")

    for c in 1:copCF.m
        obj = string(obj, "\ng cell", c, "\n")
        for f in copCF[c, :].nzind
            triangles = triangulated_faces[f]
            for tri in triangles
                t = copCF[c, f] > 0 ? tri : tri[end:-1:1]
                obj = string(obj, "f ", t[1], " ", t[2], " ", t[3], "\n")
            end
        end
    end

    return obj
end


"""
    obj2lar(path)

Read OBJ file at `path` and create a 2-skeleton as `Tuple{Points, ChainComplex}` from it.

This function does not care about eventual internal grouping inside the OBJ file.
"""
function obj2lar(path)
    vs = Array{Float64, 2}(0, 3)
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
                    v1 = parse(Int, split(elems[2], "/")[1])
                    v2 = parse(Int, split(elems[3], "/")[1])
                    v3 = parse(Int, split(elems[4], "/")[1])

                    e1 = sort([v1, v2])
                    e2 = sort([v2, v3])
                    e3 = sort([v1, v3])

                    if !(e1 in edges)
                        push!(edges, e1)
                    end
                    if !(e2 in edges)
                        push!(edges, e2)
                    end
                    if !(e3 in edges)
                        push!(edges, e3)
                    end

                    push!(faces, sort([v1, v2, v3]))
                end
            end
        end
    end

    return vs, build_cops(edges, faces)
end



function chainop2cells(copEV,copFE)
	EV = [SparseArrays.findnz(copEV[e,:])[1] for e=1:size(copEV,1)]
	FV = [collect(Set(vcat([EV[e] for e in SparseArrays.findnz(copFE[f,:])[1]]...))) 
			for f=1:size(copFE,1)]
	return EV,FV
end
