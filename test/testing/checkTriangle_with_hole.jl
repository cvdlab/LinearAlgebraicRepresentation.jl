using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

function build_copFE(FV::Lar.Cells, EV::Lar.Cells, signed=true)
	kEV = Lar.build_K(EV)
	kFV = Lar.build_K(FV)
	kFE = kFV * kEV'

	I = Int64[]; J = Int64[]; V = Int64[]
	for (i,j,v) in zip(findnz(kFE)...)
		if v == 2
			push!(I,i); push!(J,j); push!(V,1)
		end
	end
	copFE = sparse(I,J,V)
	for e=1:size(copFE,2)
		@assert length(findnz(copFE[:,e])[1])%2 == 0  "STOP: odd edeges in copFE column"
	end
	
	if signed
		cycles = Array{Array{Int64,1},1}[]
		kVE = kEV'
		for (f,verts) in enumerate(FV)
			edges = findnz(copFE[f,:])[1] # edges of face f
			numedges = length(edges) 
			@assert numedges==length(verts)  "STOP: different numbers of edges and verts"  
			
			# prepare dictionary
			v2es = []
			for v in verts
				ve = intersect(findnz(kVE[v,:])[1], edges)
				push!(v2es, (v,ve))
			end
			v2es = Dict(v2es)
			
			# walk the face and sign (-1) the counter-oriented edges
			loops = Array{Int64,1}[]
			while edges ≠ []
				global path = Int64[] 
				start = edges[1]
				push!(path, start) # set the first path edge
				v2 = EV[start][2]
				v0 = EV[start][1]
				while v2 ≠ v0 # for each cycle edge
					edge = path[end]
					nextedge = setdiff(v2es[v2],path)[1]
					invertsign = (EV[edge][2] ≠ EV[nextedge][1])
					if invertsign
						copFE[f,nextedge] = -1 # get the edge vertex # continue
						v2,_ = EV[nextedge]
					else
						_,v2 = EV[nextedge] # get the edge vertex # continue
					end
					push!(path,nextedge)
				end
				edges = setdiff(edges,path)
				@show path
				append!(loops, [path])
				@show loops
			end
			append!(cycles, [loops])
			@show cycles
		end
	end
	return copFE, cycles
end

(V, EV, FV) = ([0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.25 0.25 0.5 0.5 1.0 0.25 0.25 0.5 0.5 0.25 0.25 0.25 0.25 0.5 0.5 0.5 0.5; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.25 0.75 0.25 0.75 1.0 0.25 0.75 0.25 0.75 0.25 0.25 0.75 0.75 0.25 0.25 0.75 0.75; 0.0 1.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 1.0 -0.25 1.25 -0.25 1.25 -0.25 1.25 -0.25 1.25], Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4], [1, 5], [2, 6], [5, 6], [3, 7], [5, 7], [8, 9], [8, 10], [9, 11], [10, 11], [4, 12], [6, 12], [13, 14], [13, 15], [14, 16], [15, 16], [7, 12], [8, 17], [8, 13], [13, 18], [17, 19], [18, 20], [9, 19], [9, 14], [14, 20], [17, 21], [18, 22], [10, 21], [10, 15], [15, 22], [19, 23], [21, 23], [20, 24], [22, 24], [11, 23], [11, 16], [16, 24] ], #,    [4,14],[3,9]], 
Array{Int64,1}[[4, 2, 3, 1], [2, 5, 6, 1], [7, 9, 10, 3, 11, 5, 8, 1], [9, 10, 11, 8], [4, 13, 14, 2, 16, 15, 6, 12], [13, 14, 16, 15], [7, 4, 3, 12], [7, 5, 6, 12], [9, 19, 17, 8], [9, 13, 14, 8], [20, 13, 14, 18], [10, 17, 8, 21], [13, 10, 8, 15], [13, 15, 22, 18], [23, 19, 17, 21], [20, 24, 22, 18], [23, 9, 19, 11], [9, 14, 16, 11], [24, 14, 16, 20], [23, 10, 11, 21], [10, 16, 11, 15], [16, 15, 22, 24]])

Plasm.view(Plasm.numbering(.25)((V,[[[k] for k=1:size(V,2)],EV,FV])))

points_map = [2,4,6,12, 13,14,15,16]
edges_list = hcat([[2,6],[2,4],[6,12],[4,12],  [13,14],[13,15],[14,16],[15,16]]...)'
points =  hcat([V[:,w] for w in points_map]...)'
edges_boundary = [true,true,true,true,true,true,true,true]
points_map = vcat(points_map...)
TV = Triangle.constrained_triangulation(points,points_map,edges_list,edges_boundary)

points_map = [[2,4,6,12], [13,14,15,16]]
for triangle in TV
	for cycle in points_map[2:end]
		if setdiff(triangle, cycle)==[] # triangle ⊆ interior_cycle	
			println(triangle)
		end
	end
end
