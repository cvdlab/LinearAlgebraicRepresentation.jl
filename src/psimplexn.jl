#####################################################
### From Julia 0.6 to Julia 1.0
### --- SIMPLEXN - parallel version ---
### Authors: Fabio Fatelli, Emanuele Loprevite
#####################################################

# Generation of the output model vertices in a multiple extrusion of a LAR model
using SharedArrays

@everywhere function plarExtrude1(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}, pattern::Array{Int64,1}) where T<:Real
	V, FV = model
	d, m = length(FV[1]), length(pattern)
	dd1 = d*(d+1)
	coords = cumsum(append!([0],abs.(pattern))) # built-in function cumsum
	outVertices = @spawn [vcat(v,z) for z in coords for v in V]
	offset, outcells, rangelimit = length(V), SharedArray{Int64}(m,dd1*length(FV)), d*m
	@sync @distributed for j in 1:length(FV)
		tube = [v+k*offset for k in 0:m for v in FV[j]]
		celltube = Int64[]
		celltube = @sync @distributed (append!) for k in 1:rangelimit
			tube[k:k+d]
		end
		outcells[:,(j-1)*dd1+1:(j-1)*dd1+dd1] = reshape(celltube,dd1,m)'
	end
	p = convert(SharedArray,findall(x->x>0,pattern))
	cellGroups = SharedArray{Int64}(length(p),size(outcells)[2])
	@sync @distributed for k in 1:length(p)
		cellGroups[k,:] = outcells[p[k],:]
	end
	cellGroupsL = vec(cellGroups')
	outCellGroups = Array{Int64,1}[]
	outCellGroups = @distributed (append!) for k in 1:d+1:length(cellGroupsL)
			[cellGroupsL[k:k+d]]
		end
	return fetch(outVertices), outCellGroups
end

# Generation of simplicial grids of any dimension and shape
@everywhere function plarSimplexGrid1(shape::Array{Int64,1})
	model = [Int64[]],[[0]] # the empty simplicial model
	for item in shape # no parallel
		model = plarExtrude1(model,repeat([1],item))
	end
	return model
end

# Extraction of non-oriented (d-1)-facets of d-dimensional simplices
@everywhere using Combinatorics # for combinations() function

@everywhere function plarSimplexFacets(simplices::Array{Array{Int64,1},1})
	out = Array{Int64,1}[]
	d = length(simplices[1])
	out = @distributed (append!) for simplex in simplices
			collect(combinations(simplex,d-1))
		end
	return sort!(unique(out),lt=isless)
end

# Transformation to triangles by sorting circularly the vertices of faces
using LinearAlgebra

@everywhere function pquads2tria(model::Tuple{Array{Array{T,1},1},Array{Array{Int64,1},1}}) where T<:Real
	V, FV = model
	if typeof(V) != Array{Array{Float64,1},1}
		V = convert(Array{Array{Float64,1},1},V)
	end
	out = Array{Int64,1}[]
	nverts = length(V)-1
	for face in FV # no parallel
		arr = [V[v+1] for v in face]
		centroid = sum(arr)/length(arr)
		append!(V,[centroid])
		nverts += 1
		v1, v2 = V[face[1]+1]-centroid, V[face[2]+1]-centroid
		v3 = cross(v1,v2)
		if norm(v3) < 1/(10^3)
			v1, v2 = V[face[1]+1]-centroid, V[face[3]+1]-centroid
			v3 = cross(v1,v2)
		end
		transf = inv(copy(hcat(v1,v2,v3)'))
		verts = [(V[v+1]'*transf)'[1:end-1] for v in face]
		tcentroid = sum(verts)/length(verts)
		tverts = pmap(x->x-tcentroid,verts)
		iterator = collect(zip(tverts,face))
		rverts = [[atan(reverse(iterator[i][1])...),iterator[i][2]] for i in 1:length(iterator)]
		rvertsS = sort(rverts,lt=(x,y)->isless(x[1],y[1]))
		ord = [pair[2] for pair in rvertsS]
		append!(ord,ord[1])
		edges = [[i[2],ord[i[1]+1]] for i in enumerate(ord[1:end-1])]
		triangles = pmap(x->prepend!(x,nverts),edges)
		append!(out,triangles)
	end
	return V, out
end
