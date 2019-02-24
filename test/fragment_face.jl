using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using NearestNeighbors





#	transform sigma and related faces in space_idx
#	function face_mapping(V, FV, sigma, err=LinearAlgebraicRepresentation.ERR )
#		sigmavs = FV[sigma]; i = 1
#		compute affinely independent triple
#		while (-err < det(V[:, sigmavs[i:i+2]]) < err) && (i < length(sigmavs)-2)
#			i += 1
#		end
#		sigma translation
#		T = Lar.t(-V[:,sigmavs[i+1]]...)
#		W = T * [V[:,sigmavs] ; ones(1,length(sigmavs))]
#		t_vs = W[1:3,:]
#		translated sigma rotation
#		r1,r2,r3 = t_vs[:,i], t_vs[:,i+2], cross(t_vs[:,i], t_vs[:,i+2])
#		R = eye(4)
#		R[1:3,1:3] = inv([r1 r2 r3])
#		mapping = R * T
#		return mapping
#	end


# transform sigma and related faces in space_idx
function face_mapping(V, FV, sigma, err=LinearAlgebraicRepresentation.ERR )
@show sigma
	vs = FV[sigma]; i = 1
	# compute affinely independent triple
	n = length(vs)
	succ(i) = i % n + 1
	a = normalize(V[:,vs[succ(i)]] - V[:,vs[i]])
	b = normalize(V[:,vs[succ(succ(i))]] - V[:,vs[i]])
	c = cross(a,b)
	while (-err < det([a b c]) < err) 
		i += 1
		a = normalize(V[:,vs[succ(i)]] - V[:,vs[i]])
		b = normalize(V[:,vs[succ(succ(i))]] - V[:,vs[i]])
		c = cross(a,b)
	end
	# sigma translation
	T = Lar.t(-V[:,vs[succ(i)]]...)
	W = T * [V[:,vs] ; ones(1,length(vs))]
	t_vs = W[1:3,:]
	# translated sigma rotation
	r1,r2,r3 = t_vs[:,i], t_vs[:,succ(succ(i))], cross(t_vs[:,i], t_vs[:,succ(succ(i))])
	R = eye(4)
	R[1:3,1:3] = inv([r1 r2 r3])
	mapping = R * T
	return mapping
end


function polyline(V::Lar.Points)::Lar.LAR
	return V,[[k,k+1] for k=1:size(V,2)-1]
end

function sigmamodel(V, copEV, FV, copFE, sigma, space_idx)
	bigpi = space_idx[sigma]
	Q = Lar.face_mapping(V', FV, sigma)
	EV = [findnz(copEV[e,:])[1] for e=1:size(copEV,1)]
	FE = [findnz(copFE[e,:])[1] for e=1:size(copFE,1)]
	sigma_edges = EV[FE[sigma]]
	sigma_verts = union(sigma_edges...)
	sigma_vdict =  Dict(zip(sigma_verts, 1:length(sigma_verts)))
	sigma_lines = [[sigma_vdict[v] for v in edge] for edge in sigma_edges]
	# sigma face on z=0 plane
	S = V'[:,sigma_verts] # sigma vertices
	Z = (Q * [S; ones(1,size(S,2))])[1:3,:] # sigma mapped in z=0
	v,ev = Z,sigma_lines
	#Plasm.view(Plasm.numbering(0.5)((v,[[[k] for k=1:size(v,2)],ev])))
	return v,ev,Q
end

function sigma_intersect(V, copEV, copFE, sigma, Q, bigpi)
	EV = [findnz(copEV[e,:])[1] for e=1:size(copEV,1)]
	FE = [findnz(copFE[f,:])[1] for f=1:size(copFE,1)]
	sigma_edges = EV[FE[sigma]]
	sigma_verts = union(sigma_edges...)
	sigma_vdict =  Dict(zip(sigma_verts, 1:length(sigma_verts)))
	sigma_lines = [[sigma_vdict[v] for v in edge] for edge in sigma_edges]
	# sigma face on z=0 plane
	S = V[:,sigma_verts] # sigma vertices
	Z = (Q * [S; ones(1,size(S,2))])[1:3,:] # sigma mapped in z=0
	#Plasm.view(Z, sigma_lines)
	# initialization of line storage (to be intersected)
	linestore = [[Z[:,v1] Z[:,v2]] for (v1,v2) in sigma_lines]
	b_linenum = length(sigma_lines)
	# reindexing of ∏(σ) chain 
	bigpi_edges = [EV[fe] for fe in FE[bigpi]]
	bigpi_verts = union(union(bigpi_edges...)...)
	bigpi_vdict =  Dict(zip(bigpi_verts, 1:length(bigpi_verts)))
	bigpi_lines = [[sort([bigpi_vdict[v] for v in edge]) for edge in faceedges] 
		for faceedges in bigpi_edges]
	# bigpi trasformed in 3D according to Q mapping
	P = V[:,bigpi_verts] # bigpi vertices
	W = (Q * [P; ones(1,size(P,2))])[1:3,:] # bigpi mapped by Q
	#Plasm.view(W, union(bigpi_lines...))
	# filter on bigpi_lines that do not cross z=0
	filtered_edges = [[[v1,v2] for (v1,v2) in face_lines 
		if (sign(W[3,v1]) * sign(W[3,v2])) <= 0] 
			for face_lines in bigpi_lines]
	filtered_edges = [edge for edge in filtered_edges if edge!=[]]
	#Plasm.view(W, union(filtered_edges...))
	#Plasm.view(Plasm.numbering()((W,[[[k] for k=1:size(W,2)],union(filtered_edges...)])))
	# computation of ordered z=0 points by face
	facepoints = []
	for face in filtered_edges
		points = []
		for (v1,v2) in face
			#v1[3] + (v2[3]-v1[3])*t = 0 # compute parameter of intersection point with z=0
			if (W[3,v2] - W[3,v1]) != 0
				t = -W[3,v1] / (W[3,v2] - W[3,v1])  # z1 + (z2-z1)*t = 0
				point = W[:,v1] + (W[:,v2] - W[:,v1]) * t
				push!(points,point)
			end
		end
		if length(Set(points)) > 1
			push!(facepoints,Set(points))
		end
		# TODO: linearly order facepoints[face]
	end
	# compute z_lines
	c = collect
	z0_lines = [[[c(ps)[k] c(ps)[k+1]] for k=1:2:(length(ps)-1)] for ps in facepoints]
	
	# intersecting lines upload in linestore
	linestore = [[Z[:,v1] Z[:,v2]] for (v1,v2) in sigma_lines]
	if z0_lines != []
		linestore = append!(linestore, union(z0_lines...))
	end
	return linestore,b_linenum
end


function computeparams(linestore,b_linenum)
	m = length(linestore)
	params = [[] for i=1:m]
	for h=1:(m-1)
		(x1,y1,x2,y2) = reshape(linestore[h][1:2,:],1,4)
		for k=(h+1):m
			(x3,y3,x4,y4) = reshape(linestore[k][1:2,:],1,4)
			# intersect lines h,k
			# http://www.cs.swan.ac.uk/~cssimon/line_intersection.html
			det = (x4-x3)*(y1-y2)-(x1-x2)*(y4-y3)
			α,β = 99,99
			if det != 0.0
				a = 1/det
				b = [y1-y2 x2-x1; y3-y4 x4-x3]  # x1-x2 => x2-x1 bug in the source link !!
				c = [x1-x3; y1-y3]
				(β,α) = a * b * c
			else
				if (y1==y2) == (y3==y4) # segments collinear
					α = -x1/(x2-x1)
					β = -x3/(x4-x3)
				elseif (x1==x2) == (x3==x4)
					α = -y1/(y2-y1)
					β = -y3/(y4-y3)
				else
					 # segments parallel: no intersection
					 (β,α) = 999
				end
			end
			if 0<=α<=1 && 0<=β<=1
				push!(params[h], α)
				push!(params[k], β)
			end
		end
	end
	out = []
	for line in params
		push!(line, 0.0)
		push!(line, 1.0)
		line = sort(collect(Set(line)))
		push!(out, line)
	end	
	return out
end

function fragface(V, EV, FV, FE, space_idx, sigma)
	bigpi = space_idx[sigma]
	Q = Lar.face_mapping(V', FV, sigma)
	linestore, b_linenum = Lar.sigma_intersect(V', EV, FE, sigma, Q, bigpi)
	lineparams = Lar.computeparams(linestore,b_linenum)
	pairs = collect(zip(lineparams,linestore))

	linepoints = []
	number = 0
	for (k,(params,line)) in enumerate(pairs)
		v1,v2 = line[1:2,1],line[1:2,2]
		if params==[0.0, 1.0]
			points = [v1,v2]
			if k<=b_linenum number += 1 end
		elseif length(params)>2 
			points = [v1+α*(v2-v1) for α ∈ params]
			if k<=b_linenum number += length(params)-1 end
		end
		push!(linepoints,points)
	end
	#linepoints = collect(Set(linepoints))
	verts = []
	edges = Array{Array{Int},1}()
	offset = 0
	for (h,pointline) in enumerate(linepoints)
		append!(verts, pointline)
		append!(edges, [[k,k+1] for k=1:length(pointline)-1]+offset)
		offset = length(verts)
	end
	verts = hcat(verts...)
	edges = convert(Lar.Cells, edges)
	return verts, edges, number
end


function removevertices(verts, edges, err=LinearAlgebraicRepresentation.ERR)
	vertices = collect(Set(vcat(edges...)))
	vertexnum = length(vertices)
	vdict = Dict(zip(vertices,1:vertexnum))
	cells0D = cells0D = zeros(2,vertexnum)
	cells1D = Array{Int64,1}[]
	for (v1,v2) in edges
		push!(cells1D,[vdict[v1],vdict[v2]])
	end
	#vdict = Dict(zip(values(vdict),keys(vdict))) # inverted vdict
	for key in keys(vdict)
		cells0D[:,vdict[key]] = verts[:,key]
	end
	return cells0D, cells1D
end
	
function merge_vertices(rV::Lar.Points, rEV::Lar.ChainOp, rFE::Lar.ChainOp, 
err=Lar.ERR)
	verts = rV'
	edges = [findnz(rEV[e,:])[1] for e=1:size(rEV,1)]
	faces = [findnz(rFE[f,:])[1] for f=1:size(rFE,1)]


	function make_quotients(verts, faces, edges, err=Lar.ERR)
print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
@show verts
@show edges
@show faces
PRECISION = 5.0

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
			key = map(Lar.approxVal(PRECISION), nverts[:,k])
			vertdict[key] = k
		end

		nedges = Array{Array{Int},1}()
		for (v1,v2) in edges
			(w1,w2) = abs.((v1,v2))
			key1 = map(Lar.approxVal(PRECISION), verts[:,w1])
			key2 = map(Lar.approxVal(PRECISION), verts[:,w2])
			push!(nedges, [vertdict[key1], vertdict[key2]])
		end
		nedges = Set(map(sort,nedges))
		nedges = sort(collect(nedges), by=x->(x[1],x[2]))
		nedges = filter( x->x[1]!=x[2], nedges )	
#		nedges = sort(nedges, by=x->(x[1],x[2]))
#		nedges = sort(collect(Set(nedges)), by=x->(x[1],x[2]))
#		nedges = filter( x->x[1]!=x[2], nedges )
#		nedges = convert(Array{Array{Int64,1},1}, nedges)	

		nfaces = Array{Array{Int},1}()
		for (k,face) in enumerate(faces)
			nface = Int64[]
			for edge in face
				(v1,v2) = edges[edge]
				(w1,w2) = abs.((v1,v2))
				key1 = map(Lar.approxVal(PRECISION), verts[:,w1])
				key2 = map(Lar.approxVal(PRECISION), verts[:,w2])
				push!(nface, vertdict[key1])
				push!(nface, vertdict[key2])
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

print("before make_quotients ^^^^^^^^^^^^^^^^^^^")
	V,FV,EV = make_quotients(verts, faces, edges, err)
	return V,FV,EV	
end
	
function (verts, edges, err=LinearAlgebraicRepresentation.ERR)
	vertsnum = size(verts, 2)
	kdtree = KDTree(verts)
	todelete = []
	newverts = zeros(Int, vertsnum)
	i = 1
	for vi in 1:vertsnum
	   if !(vi in todelete)
		   nearvs = inrange(kdtree, verts[:,vi], err)
		   newverts[nearvs] = i
		   nearvs = setdiff(nearvs, vi)
		   todelete = union(todelete, nearvs)
		   i = i + 1
	   end
	end
	nverts = verts[:,setdiff(collect(1:vertsnum), todelete)]
	vertdict = DataStructures.OrderedDict()
	for k=1:size(nverts,2)
		key = map(Lar.approxVal(PRECISION), nverts[:,k])
		vertdict[key] = k
	end
	nedges = Array{Array{Int},1}()
	for (v1,v2) in edges
		(w1,w2) = abs.((v1,v2))
		key1 = map(Lar.approxVal(PRECISION), verts[:,w1])
		key2 = map(Lar.approxVal(PRECISION), verts[:,w2])
		push!(nedges, [vertdict[key1], vertdict[key2]])
	end
	nedges = Set(map(sort,nedges))
	nedges = sort(collect(nedges), by=x->(x[1],x[2]))
	nedges = filter( x->x[1]!=x[2], nedges )	
	return nverts,nedges
end

function preprocessing(V,copEV,FV)
	EV = [findnz(copEV[e,:])[1] for e=1:size(copEV,1)]
	#copEV = Lar.coboundary_0(EV)
	#copFE = build_copFE(FV, EV)
	copFE = Lar.build_copFE(FV, EV)
	copEF = copFE'
	EF = [findnz(copEF[e,:])[1] for e=1:size(copEF,1)]
	FE = [findnz(copFE[e,:])[1] for e=1:size(copFE,1)]

	sp_idx = Lar.Arrangement.spatial_index(V, copEV, copFE)
	space_idx = Array{Int,1}[]
	for sigma = 1:length(FV)
		edges = FE[sigma];
		# remove faces incident on sigma edges
		facesofedges = Set(vcat([EF[edge] for edge in edges]...));
		push!(space_idx, setdiff(sp_idx[sigma], facesofedges))
	end
	return V, copEV, copFE, space_idx
end

function computefragments(V, EV, FV, FE, space_idx, sigma)
	#v,ev = Lar.sigmamodel(V, EV, FE, sigma, space_idx)
	#Plasm.view(Plasm.numbering(0.4)((v,[ [[k] for k=1:size(v,2)],ev ])))

	verts, edges, b_linenum = Lar.fragface(V, EV, FV, FE, space_idx, sigma) #OK
	v,ev,Q = Lar.sigmamodel(V, EV, FV, FE, sigma, space_idx) #OK
	classify = Lar.pointInPolygonClassification(v,ev)
	internal = []
	for (k,edge)  in enumerate(edges)
		v1,v2 = edge
		if k <= b_linenum
			push!(internal, edge)
		else
			vm = (verts[:,v1]+verts[:,v2])/2
			if (classify(vm)=="p_in")
				push!(internal, edge)
			end
		end
	end
	
	nverts, nedges = Lar.removevertices(verts, internal)
	w, ev = Lar.mergevertices(nverts, nedges)
	cop_ev = Lar.build_copEV(ev,true)
	cop_fe = Lar.Arrangement.minimal_2cycles(w',cop_ev)	#OK
	shell_num = Lar.Arrangement.get_external_cycle(w', cop_ev, cop_fe)
	for j=1:cop_fe.n  
		cop_fe[shell_num,j]=0 
	end
	cop_fe = dropzeros(cop_fe)
	faces = [findnz(cop_fe[e,:])[1] for e=1:size(cop_fe,1) if e≠shell_num]
	fv = [collect(Set(vcat([ev[e] for e in face]...))) for face in faces]
	vv = [[k] for k=1:size(w,2)]
	#Plasm.view(Plasm.numbering(0.4)((w, [vv, ev, fv])))
		
	w = [ w; zeros(size(w,2))'; ones(size(w,2))' ]
	v = (inv(Q) * w)[1:3,:]
	vv = [[k] for k=1:size(w,2)]
	#Plasm.view(Plasm.numbering(0.4)((v, [vv, ev, fv])))
	return v', Lar.coboundary_0(ev), cop_fe
end









