# https://www.cp.eng.chula.ac.th/~attawith/class/linint
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, DataStructures


function cuboids(n,scale=1.; grid=[1,1])
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		dims = rand(Float64, 2)
		if grid==[1,1]
			V,(_,EV,_) = Lar.cuboid(corner,true,corner+dims)
		else
			V,(_,EV,_) = Lar.cuboidGrid(grid, true)
			g = Lar.Struct([ Lar.t(corner...), Lar.s(dims...), (V,EV) ])
			V,EV = Lar.struct2lar(g)
		end
		center = (corner + corner+dims)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end

function presorted(V,EV)
	lines = [[V[:,u],V[:,v]] for (u,v) in EV]
	orientedlines = Array{Array{Float64,1},1}[]
	for (v1,v2) in lines
		if v1 > v2 
			v1,v2 = v2,v1
		end
		push!(orientedlines, [v1,v2])
	end
	outlines = sort(orientedlines)
end

function sweepline(V,EV)
	segments = presorted(V,EV)
	evpairs = [[(v1,"start",k), (v2,"end",k)] for (k,(v1,v2)) in enumerate(segments)]
	events = sort(cat(evpairs))
	# Initialize event queue ξ = all segment endpoints Sort ξ by increasing x and y
	pqvalues = events
	pqkeys = [(e[1],e[3]) for e in events]
	ξ = PriorityQueue(zip(pqkeys,pqvalues))
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Int,Array{Float64,1}}()	
	i = startof(SL)
	# Initialized output intersection list Λ to be empty
	Λ = Array{Int64,1}[]
	while length(ξ) ≠ 0 # (ξ is nonempty)
		E = peek(ξ).second # the next event (k => v) from ξ 
		dequeue!(ξ) 
		if E[2] == "start"# (E is a left endpoint)
			# segE = E's segment
			e = E[3]; segE = EV[e]
			# Add segE to SL
			i = insert!(SL, e, segE)
			if segA ≠ beforestartsemitoken(SL)
				# segA = the segment above segE in SL 
				segA = regress((SL,i))
				# (I = Intersect( segE with segA) exists)
				# Insert I into ξ
			end
			if deref(SL,i) ≠ last(SL)
				# segB = the segment below segE in SL 
				segB = advance((SL,i))
				# (I = Intersect( segE with segB) exists)
				# Insert I into ξ
			end
		elseif E[2] == "end" # (E is a right endpoint)
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


# EXAMPLE 1
# data generation
V,EV = cuboids(10, .5)
V = Plasm.normalize(V,flag=true)
model2d = V,EV
Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))



#===================================================================
#===================================================================
#===================================================================
#===================================================================



# https://www.cp.eng.chula.ac.th/~attawith/class/linint
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, DataStructures


function cuboids(n,scale=1.; grid=[1,1])
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		dims = rand(Float64, 2)
		if grid==[1,1]
			V,(_,EV,_) = Lar.cuboid(corner,true,corner+dims)
		else
			V,(_,EV,_) = Lar.cuboidGrid(grid, true)
			g = Lar.Struct([ Lar.t(corner...), Lar.s(dims...), (V,EV) ])
			V,EV = Lar.struct2lar(g)
		end
		center = (corner + corner+dims)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end


function presorted(V,EV, preview=false)
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
	W,EW  = lines2lar(newlines)
	if preview
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



"""
	intersection(V::Lar.Points, EV::Lar.Cells)
		(e1::Int, e2::Int)::Union{Nothing, Array}
		
Compute the intersection point of two line segments `e1` and `e2`,
of extreme point indices `EV[e1]` and `EV[e2]`, with coordinates 
in `V` array.
"""
function intersection(V::Lar.Points, EV::Lar.Cells)
	function intersection0(e1::Int, e2::Int)::Union{Nothing, Array}
		x1,y1,x2,y2 = reshape(V[:,EV[e1]], 4,1)
		x3,y3,x4,y4 = reshape(V[:,EV[e2]], 4,1)
	
		# intersect lines e1,e2
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
			p = v1 + α * (v2 - v1)
			return p
		end
	end
	return intersection0
end


#=
function sweepline(V,EV) # input: unordered LAR
	# data preprocessing
	V,EV = presorted(V,EV) # V and EV forward ordered in x,y 
	vertdict = Dict([(v,i) for (i,v) in enumerate([V[:,k] for k=1:size(V,2)])])
	startVE = Array{Any,1}[]; for k=1:size(V,2) push!(startVE,[]) end
	endVE = Array{Any,1}[]; for k=1:size(V,2) push!(endVE,[]) end
	for (e,(u,v)) in enumerate(EV) 
		push!(startVE[u],e)
		push!(endVE[v],e)
	end
	# Initialize event priority queue ξ = all segment endpoints; 
	pqvalues = [V[:,k] for k=1:size(V,2)]
	pqkeys = [x for (x,y) in pqvalues]
	# ξ sorted by increasing x and y
	## pqkeys: used float values in order to orderly insert in ξ
	## x-values of discovered intersection points 
	ξ = PriorityQueue(zip(pqkeys,pqvalues))
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Float64,Int}()	
	i = startof(SL)
	# Initialized output intersection list Λ to be empty
	Λ = Array{Int64,1}[]
	while length(ξ) ≠ 0 # (ξ is nonempty)
	
		E = peek(ξ) # the next event (k => v) from ξ 
		dequeue!(ξ)
		v = vertdict[E[2]] # event node
		
		for e in union(startVE[v], endVE[v]) # edges incident on event node 
			# Add segE to SL
			i = insert!(SL, v, segE)
				
			if e in startVE[v] # (E is a left endpoint)
				segE = e # segE = E segment
				# segA = the segment above segE in SL 
				segA = regress((SL,i))
				# segB = the segment below segE in SL 
				segB = advance((SL,i))
				if above ≠ afterendsemitoken 
					# if I = Intersect( segE with segA ) exists
					I = intersection(V,EV)(segE,segA)
					if I ≠ None
						segE1,segE2 = segE,segA
						# Insert I into ξ
						key = I[1]
						val = I
						enqueue!(ξ, key, val) 
					end
				end
				I = intersection(V,EV)(segE,segB)
				# (I = Intersect( segE with segB) exists)
				if I ≠ None
					segE1,segE2 = segE,segB
					# Insert I into ξ
					key = I[1]
					val = I
					enqueue!(ξ, key, val) 
				end
				
			elseif e in endVE[v] # (E is a right endpoint)  ## KO => may be both !!!
				segE = e # segE = E segment
				# segA = the segment above segE in SL 
				segA = regress((SL,i))
				# segB = the segment below segE in SL 
				segB = advance((SL,i))
				if above ≠ afterendsemitoken 
					# if I = Intersect( segE with segB ) exists
					I = intersection(V,EV)(segE,segB)
					if I ≠ None
						segE1,segE2 = segE,segB
						# Insert I into ξ
						key = I[1]
						val = I
						enqueue!(ξ, key, val) 
					end
				end
				
			else # E is an intersection event
				@assert E[2] !∈ keys(vertdict) # carefull: bad complexity
				push!(Λ, E[2]) # Add E to the output list Λ of intersection points
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
		end
		# remove E from ξ
	end # ξ is empty
	return Λ
end
=#

#=
# EXAMPLE 1
# data generation
V,EV = cuboids(10, .5)
V = Plasm.normalize(V,flag=true)
model2d = V,EV
Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))
# data sorting 
W,EW = presorted(V,EV)
Plasm.view(Plasm.numbering(.15)((W,[[[k] for k=1:size(W,2)], EW ])))


# EXAMPLE 2
# data generation
V,(_,EV,FV) = Lar.cuboidGrid([3,3],true)
Plasm.view(Plasm.numbering(1)((V,[[[k] for k=1:size(V,2)], EV ])))
# data sorting 
W,EW = presorted(V,EV)
Plasm.view(Plasm.numbering(1)((W,[[[k] for k=1:size(W,2)], EW ])))
=#


# EXAMPLE 3
# data generation
V,EV = cuboids(4, .5; grid=[2,2])
V = Plasm.normalize(V,flag=true)
Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV ])))
# data sorting 
W,EW = presorted(V,EV)
Plasm.view(Plasm.numbering(.15)((W,[[[k] for k=1:size(W,2)], EW ])))



    # Test all methods of SortedMultiDict except loops


julia> smd = SortedMultiDict(o)
SortedMultiDict(Base.Order.ForwardOrdering(),)

julia> push!(smd, 3.5=>35)
SortedMultiDict(Base.Order.ForwardOrdering(),3.5 => 35)

julia> push!(smd, 1.5=>15)
SortedMultiDict(Base.Order.ForwardOrdering(),1.5 => 15, 3.5 => 35)

julia> push!(smd, 3.5=>350)
SortedMultiDict(Base.Order.ForwardOrdering(),1.5 => 15, 3.5 => 35, 3.5 => 350)

julia> first(smd)
1.5 => 15

julia> last(smd)
3.5 => 350

SortedMultiDict(o, [(1.5,15),(2.5,25),(3.5,35)])


compare(sc, st1, st2)

advance((sc,st))


regress((sc,st))



    factors = SortedMultiDict{Float64,Int}()
    N = 1000
    len = 0
    sum1 = 0
    sum2 = 0
    for factor = 1 : N
        for multiple = factor : factor : N
            insert!(factors, multiple, factor)
            global sum1 += multiple
            global sum2 += factor
            global len += 1
        end
    end



for st in eachindex(inclusive(SL,searchsortedfirst(SL, 1.),
                                searchsortedfirst(SL, 2.)))
@show st  end


startof(SL)
