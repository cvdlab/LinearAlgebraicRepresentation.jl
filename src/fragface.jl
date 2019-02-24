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

"""

Intersect `sigma` edges with edges in `bigpi`.

# Example 2D

```julia


```
"""
function sigma_intersect(V, EV, FV, sigma, Q, bigpi)

	cscEV = Lar.characteristic_matrix(EV)

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



"""

Intersect pairs of `line` in `linestore`.

"""
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


"""


"""
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

