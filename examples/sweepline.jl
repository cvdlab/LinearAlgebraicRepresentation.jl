# https://www.cp.eng.chula.ac.th/~attawith/class/linint
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, DataStructures

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
		x1,y1,x2,y2 = cat(lines[h])
		
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

function presortedlines(V,EV)
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


"""
	intersection(V::Lar.Points, EV::Lar.Cells)
		(e1::Int, e2::Int)::Union{Nothing, Array}
		
Compute the intersection point of two line segments `line1` and `line2`,
of extreme point indices `EV[e1]` and `EV[e2]`, with given coordinates 
in `V` array.

# Example

```julia
julia> line1 = V[:,EV[14]]
2×2 Array{Float64,2}:
 0.0379479  0.168838
 0.90315    1.0     

julia> line2 = V[:,EV[23]]
2×2 Array{Float64,2}:
 0.0307485  0.312532
 0.93606    0.957881

julia> intersection(line1,line2)
2-element Array{Float64,1}:
 0.08846469268137241
 0.9405292234774614 

julia> intersection(line1,line2)
2-element Array{Float64,1}:
 0.08846469268137241
 0.9405292234774614 

julia> intersection(line1,-1*line2) # returns a `Nothing` value

julia> 
```
"""
function intersection(line1,line2)::Union{Nothing, Array}
	x1,y1,x2,y2 = reshape(line1, 4,1)
	x3,y3,x4,y4 = reshape(line2, 4,1)

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
		p = [x1,y1] + α * ([x2,y2] - [x1,y1])
		return p
	end
end


#=
function sweepline(V,EV)
	segments = presortedlines(V,EV)
	evpairs = [[(v1,"start",k), (v2,"end",k)] for (k,(v1,v2)) in enumerate(segments)]
	events = sort(cat(evpairs))
	# Initialize event queue ξ = all segment endpoints Sort ξ by increasing x and y
	pqvalues = events
	pqkeys = [(e[1],e[3]) for e in events]
	ξ = PriorityQueue(zip(pqkeys,pqvalues))
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Int,Array{Float64,1}}()
	@assert isempty(SL) # debug	
	@assert beforestartsemitoken(SL)==DataStructures.Tokens.IntSemiToken(1) # debug
	@assert pastendsemitoken(SL)==DataStructures.Tokens.IntSemiToken(2) # debug
	
	i = startof(SL)
	@assert i = startof(SL)==DataStructures.Tokens.IntSemiToken(2) # debug
	# Initialized output intersection list Λ to be empty
	Λ = Array{Int64,1}[]
	while length(ξ) ≠ 0 # (ξ is nonempty)
		E = peek(ξ).second; dequeue!(ξ)  # the next v event (k => v) from ξ 
		
		if E[2] == "start"# (E is a left endpoint)
			# segE = E's segment
			e = E[3]; segE = EV[e] # ==> segE = V[:,EV[e]] is better ??
			# Add segE to SL
			i = insert!(SL, e, segE)  #  SortedMultiDict, key, value
			if first(SL) !== last(SL)
				# segA = the segment above segE in SL 
				j = regress((SL,i))
				E = peek(ξ).second; a = E[3]; segA = EV[a]
				I = intersection(segE,segA)
				global seg1 = segE; global seg2 = segA
				# (if Intersect( segE with segA) exists)
				if typeof(I) ≠ Nothing
					# Insert I into ξ
					key = (E[1],E[3]); val = (I,"int",(e,a))
					enqueue!(ξ, key, val) 
				end
			end
			if i < pastendsemitoken(SL)
				# segB = the segment below segE in SL 
				k = advance((SL,i))
				E = peek(ξ).second; b = E[3]; segB = EV[b]
				I = intersection(segE,segB)
				global seg1 = segE; global seg2 = segB
				# (if Intersect( segE with segB) exists)
				if typeof(I) ≠ Nothing
					# Insert I into ξ
					key = (E[1],E[3]); val = (I,"int",(e,b))
					enqueue!(ξ, key, val) 
			end
		elseif E[2] == "end" # (E is a right endpoint)
			# segE = E's segment
			e = E[3]; segE = EV[e]
			# segA = the segment above segE in SL 
			if segE ≠ beforestartsemitoken(SL)
				# segA = the segment above segE in SL 
				j = regress((SL,i))
				E = peek(ξ).second; a = E[3]; segA = EV[a]
			end
			# segB = the segment below segE in SL 
			if deref((SL,i)) ≠ last(SL)
				# segB = the segment below segE in SL 
				k = advance((SL,i))
				E = peek(ξ).second; b = E[3]; segB = EV[b]
			end
			# Remove segE from SL
			delete!((SL,i))
			# (I = Intersect( segA with segB) exists)
			I = intersection(segA,segB)
			global seg1 = segA; global seg2 = segB
			if  typeof(I) ≠ Nothing
				h = searchsortedfirst(LS, k)
				# (I is not in ξ already) 
				if h = pastendsemitoken(SL)
					# Insert I into ξ
					E = (I,"int",(a,b))
					key = (E[1],E[3]); val = (I,"int",(a,b))
					enqueue!(ξ, key, val) 
				end
			end
		else # E is an intersection event
			@assert E[2] == "int"
			# Add E to the output list Λ
			push!(Λ, E)
			# segE1 above segE2 be E's intersecting segments in SL 
			# segA = the segment above segE2 in SL
			a = searchsortedfirst(LS, seg2)  # TODO: check seg1's type
			# segB = the segment below segE1 in SL
			b = searchsortedlast(LS, seg1)  # TODO: check seg2's type
			# Swap their positions so that segE2 is now above segE1 
			??? HOW TO DO ???
			# (I = Intersect(segE2 with segA) exists)
			I = intersection(segE2,segA)
			if typeof(I) ≠ Nothing
				k = searchsortedfirst(LS, a)
				# (I is not in ξ already) 
				if k = beforestartsemitoken(SL) 
					# Insert I into ξ
					E = (I,"int",(k,a))
					key = (E[1],E[3]); val = (I,"int",(k,a))
					enqueue!(ξ, key, val) 
				end
			end
			I = intersection(segE1,segB)
			# (I = Intersect(segE1 with segB) exists)
			if typeof(I) ≠ Nothing
				h = searchsortedfirst(LS, b)
				# (I is not in ξ already) 
				if h = pastendsemitoken(SL) 
					# Insert I into ξ
					E = (I,"int",(h,b))
					key = (E[1],E[3]); val = (I,"int",(h,b))
					enqueue!(ξ, key, val) 
				end
			end
		end
		# remove E from ξ
		delete!(ξ, E[1])  
	end
	return Λ
end
=#

# EXAMPLE 0
# data generation
lines = [[[1,3],[10,5]],[[2,6],[11,2.5]],
[[3,4],[8,2.25]],[[5,1.25],[12,5]]]
V,EV = lines2lar(lines)
Plasm.view(Plasm.numbering(2.)((V,[[[k] for k=1:size(V,2)], EV])))
V = Plasm.normalize(V,flag=true)
W,EW = presorted(V,EV)
Plasm.view(Plasm.numbering(.25)((W,[[[k] for k=1:size(W,2)], EW ])))


# EXAMPLE 1
# data generation
V,EV = cuboids(10, .5)
V = Plasm.normalize(V,flag=true)
model2d = V,EV
Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))


	i1,i2 = searchequalrange(factors, 80)
	compare(factors,regress((factors,i)),i2) == 0
    my_assert(deref((factors,i)) == Pair(2,1))
    my_assert(deref_key((factors,i)) == 2)
    my_assert(deref_value((factors,i)) == 1)
	startof(SL)==lastindex(SL)	
    my_assert(haskey(factors, 60))
    my_assert(!haskey(factors, -1))

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
