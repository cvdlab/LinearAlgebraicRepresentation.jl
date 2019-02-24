# https://www.cp.eng.chula.ac.th/~attawith/class/linint
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, DataStructures


function cuboids(n,scale=1.)
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		dims = rand(Float64, 2)
		V,(_,EV,_) = Lar.cuboid(corner,true,corner+dims)
		center = (corner + corner+dims)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end


function linesort(V,EV, preview=true)
	lines = [[V[:,u],V[:,v]] for (u,v) in EV]
	orientedlines = Array{Array{Float64,1},1}[]
	for (v1,v2) in lines
		if v1 > v2 
			v1,v2 = v2,v1
		end
		push!(orientedlines, [v1,v2])
	end
	inputlines = sort(orientedlines)
	newlines = map(cat,inputlines)
	if preview
		W,EW  = lines2lar(newlines)
		Plasm.view(Plasm.numbering(.1)((W,[[[k] for k=1:size(W,2)], EW ])))
	end
	vsorted = sort([[W[:,k],k] for k=1:size(W,2)])
	V = hcat([vsorted[k][1] for k=1:length(vsorted)]...)
	edgedict = OrderedDict([(oldv,newv) for (newv,(vert,oldv)) in enumerate(vsorted)])
	EV = [[edgedict[v1],edgedict[v2]] for (v1,v2) in EW]
	return V,EV
end

function lines2lar(lines)
	vertdict = OrderedDict{Array{Float64,1}, Int64}()
	EV = Array{Int64,1}[]
	idx = 0
	for h=1:length(lines)
		x1,y1,x2,y2 = lines[h]
		
		if ! haskey(vertdict, [x1,y1])
			idx += 1
			vertdict[[x1,y1]] = idx
		end
		if ! haskey(vertdict, [x2,y2])
			idx += 1
			vertdict[[x2,y2]] = idx
		end
		v1,v2 = vertdict[[x1,y1]],vertdict[[x2,y2]]
		push!(EV, [v1,v2])
	end
	V = hcat(collect(keys(vertdict))...) 
	return V,EV
end

#=
function sweepline(V,EV)
	# data preprocessing
	lines = linesort(V,EV, preview=true)
	# Initialize event queue ξ = all segment endpoints; Sort ξ by increasing x and y
	# Initialize sweep line SL to be empty
	# Initialize output intersection list Λ to be empty
	while # (ξ is nonempty)
		# E = the next event from ξ 
		if # (E is a left endpoint)
			# segE = E's segment
			# Add segE to SL
			# segA = the segment above segE in SL 
			# segB = the segment below segE in SL 
			if # (I = Intersect( segE with segA) exists)
				# Insert I into ξ
			end
			if # (I = Intersect( segE with segB) exists)
				# Insert I into ξ
			end
		elseif # (E is a right endpoint)
			# segE = E's segment
			# segA = the segment above segE in SL 
			# segB = the segment below segE in SL 
			# Remove segE from SL
			if # (I = Intersect( segA with segB) exists)
				if # (I is not in ξ already) 
					# Insert I into ξ
				end
			end
		else # E is an intersection event
			# Add E to the output list Λ
			# segE1 above segE2 be E's intersecting segments in SL 
			# Swap their positions so that segE2 is now above segE1 
			# segA = the segment above segE2 in SL
			# segB = the segment below segE1 in SL
			if # (I = Intersect(segE2 with segA) exists)
				if # (I is not in ξ already) 
					# Insert I into ξ
				end
			end
			if # (I = Intersect(segE1 with segB) exists)
				if # (I is not in ξ already) 
					# Insert I into ξ
				end
			end
		end
		# remove E from ξ
	end
	return Λ
end
=#

# EXAMPLE
# data generation
V,EV = cuboids(10, .5)
V = Plasm.normalize(V,flag=true)
model2d = V,EV
Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))
# data sorting 
W,EW = linesort(V,EV)
Plasm.view(Plasm.numbering(.15)((W,[[[k] for k=1:size(W,2)], EW ])))


# EXAMPLE
# data generation
V,(_,EV,FV) = Lar.cuboidGrid([3,3],true)
Plasm.view(Plasm.numbering(1)((V,[[[k] for k=1:size(V,2)], EV ])))
# data sorting 
W,EW = linesort(V,EV)
Plasm.view(Plasm.numbering(1)((W,[[[k] for k=1:size(W,2)], EW ])))


