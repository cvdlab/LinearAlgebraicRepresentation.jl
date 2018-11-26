using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

"""
#Input:[∂d−1] # Compressed Sparse Column (CSC) signed matrix (aij), where aij ∈ {−1,0,1} 
#Output: [∂d+] # CSC signed matrix
#
#[∂d+] = [] ; m,n = [∂d−1].shape ; marks = Zeros(n) # initializations 
#while Sum(marks) < 2n do
#	σ = Choose(marks) # select the (d − 1)-cell seed of the column extraction 
#	if marks[σ] == 0 then [cd−1] = [σ]
#	else if marks[σ] == 1 then [cd−1] = [−σ]
#	[cd−2] = [∂d−1] [cd−1] # compute boundary cd−2 of seed cell 
#	while [cd−2] 􏰁≠ [] do # loop until boundary becomes empty
#		corolla = []
#		for τ ∈ cd−2 do # for each “hinge” τ cell
#			[bd−1] = [τ]^t [∂d−1] #compute the τ coboundary
#			pivot = {|bd−1|} ∩ {|cd−1|} # compute the τ support
#			if τ > 0 then adj = Next(pivot,Ord(bd−1)) # compute the new adj cell 
#			else if τ < 0 then adj = Prev(pivot,Ord(bd−1))
#			if ∂d−1[τ,adj] ≠􏰁 ∂d−1[τ,pivot] then corolla[adj] = cd−1[pivot] # orient adj
#			else corolla[adj] = −(cd−1[pivot])
#		end
#		[cd−1] += corolla # insert corolla cells in current cd−1
#		[cd−2] = [∂d−1][cd−1] # compute again the boundary of cd−1 
#	end
#	for σ ∈ cd−1 do marks[σ] += 1 # update the counters of used cells
#	[∂d+] += [cd−1] # append a new column to [∂d+] 
#end
#return [∂d+]
"""
(V, EV, FV) = ([0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.25 0.25 0.5 0.5 1.0 0.25 0.25 0.5 0.5 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.5; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.25 0.75 0.25 0.75 1.0 0.25 0.75 0.25 0.75 0.25 0.25 0.75 0.75 0.25 0.25 0.75 0.75; 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 1.0 -0.25 1.25 -0.25 1.25 -0.25 1.25 -0.25 1.25], Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4], [1, 5], [2, 6], [5, 6], [3, 7], [5, 7], [8, 9], [8, 10], [9, 11], [10, 11], [4, 12], [6, 12], [13, 14], [13, 15], [14, 16], [15, 16], [7, 12], [8, 17], [8, 13], [13, 18], [17, 19], [18, 20], [9, 19], [9, 14], [14, 20], [17, 21], [18, 22], [10, 21], [10, 15], [15, 22], [19, 23], [21, 23], [20, 24], [22, 24], [11, 23], [11, 16], [16, 24]], #,    [4,14],[3,9]], 
Array{Int64,1}[[4, 2, 3, 1], [2, 5, 6, 1], [7, 9, 10, 3, 11, 5, 8, 1], [9, 10, 11, 8], [4, 13, 14, 2, 16, 15, 6, 12], [13, 14, 16, 15], [7, 4, 3, 12], [7, 5, 6, 12], [9, 19, 17, 8], [9, 13, 14, 8], [20, 13, 14, 18], [10, 17, 8, 21], [13, 10, 8, 15], [13, 15, 22, 18], [23, 19, 17, 21], [20, 24, 22, 18], [23, 9, 19, 11], [9, 14, 16, 11], [24, 14, 16, 20], [23, 10, 11, 21], [10, 16, 11, 15], [16, 15, 22, 24]])

Plasm.view(Plasm.numbering(.25)((V,[[[k] for k=1:size(V,2)],EV,FV])))

copFE = Lar.build_copFE(V,FV,EV)

   


function area(v1,v2,v3)
	u = V[:,v2]-V[:,v1]
	v = V[:,v3]-V[:,v1]
	norm(cross(u,v)) # actually, to be divided by two
end

function interior_to_f(triangle,f,FE)
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
	checkpoint = 0.9 .* p2D[:,1] + 0.05 .* p2D[:,2] + 0.05 .* p2D[:,3]
	
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


function ordering(triangles)
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

function ord(hinge,bd1,V,FV,EV,FE)
	cells = findnz(bd1)[1]
	triangles = []
	for f in cells
		v1,v2 = EV[hinge]		
		index = findfirst(v -> (area(v1,v2,v)≠0), FV[f])
		v3 = FV[f][index]
		
		# test if [v1,v2,v3] interior to f
		while true
			if interior_to_f([v1, v2, v3],f,FE)
				push!(triangles, [v1,v2,v3])
				break
			else
				index = findnext(v -> (area(v1,v2,v)≠0), FV[f], index+1)
				v3 = FV[f][index]
			end
		end
	end
	order = ordering(triangles)
	return [cells[index] for index in order]
end


function next(cycle, pivot)
	len = length(cycle)
	ind = find(x -> x==pivot, cycle)[1]
	nextIndex = ind==len ? 1 : ind+1
	return cycle[nextIndex][1]
end

function prev(cycle, pivot)
	len = length(cycle)
	ind = find(x->x==pivot, cycle)[1]
	nextIndex = ind==1 ? len : ind-1
	return cycle[nextIndex][1]
end

function algorithm_one(V,FV,EV,copFE)
	copEF = copFE'
	FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
	# Initializations
	m,n = size(copEF)
	marks = zeros(Int,n);
	I = Int64[]; J = Int64[]; V = Int8[]; 
	jcol = 0
	choose(marks) = findfirst(x -> x<2, marks)

	# Main loop (adding one copFC's column stepwise)
	while sum(marks) < 2n
		# select a (d−1)-cell, "seed" of the column extraction
		σ = choose(marks)
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
			@show corolla
			# for each “hinge” τ cell
			for τ ∈ (.*)(findnz(cd2)...)
				@show τ
				#compute the  coboundary
				tau = sparsevec([abs(τ)], [sign(τ)], m)
				bd1 = transpose(transpose(tau) * copEF)
				@show bd1
				cells2D = findnz(bd1)[1]
				# compute the  support
				pivot = intersect(cells2D, findnz(cd1)[1])[1] # check for sign?
				# compute the new adj cell
				fan = ord(abs(τ),bd1,V,FV,EV,FE) # ord(pivot,bd1)
				if τ > 0 
					adj = next(fan,pivot)
				elseif τ < 0 
					adj = prev(fan,pivot)
				end
				@show adj
				# orient adj
				if copEF[abs(τ),adj] ≠ copEF[abs(τ),pivot] 
					corolla[adj] = cd1[pivot]
				else
					corolla[adj] = -(cd1[pivot])
				end
			end
			# insert corolla cells in current cd1
			for (k,val) in zip(findnz(corolla)...) 
				cd1[k] = val 
			end 
			@show cd1
			# compute again the boundary of cd1
			cd2 = copEF * cd1
			@show cd2
		end
		for σ ∈ findnz(cd1)[1]
			# update the counters of used cells
			marks[σ] += 1
		end
		# append a new column to [∂d+]
		# copFC += cd1
		rows, vals = findnz(cd1)
		jcol += 1
		append!(I,rows)
		append!(J,[ jcol for k=1:nnz(cd1) ])
		append!(V,vals)
	end
	copFC = sparse(I,J,V)
	return copFC
end

copFC = algorithm_one(V,FV,EV,copFE)
Matrix(copFC)



