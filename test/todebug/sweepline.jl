# https://www.cp.eng.chula.ac.th/~attawith/class/linint
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm, DataStructures


"""
	presorted(V, EV [,preview=false] )

Generate the presorting needed by the scan-line algorithm (Bentley-Ottmann, 1979).

When two or more 2D (left to right oriented) lines have the same start point, 
they are ordered according to forward values of the y-coord of the end line point.
This guarantess that each other line having the same first vertex is inserted in the right 
scan-line data structure interval, i.e. between the two lines closely above and below it.

# Example

```
julia> lines = [[Float64[1,2],Float64[4,3]],
       [Float64[1,2],Float64[4,2.999]]]
2-element Array{Array{Array{Float64,1},1},1}:
 [[1.0, 2.0], [3.0, 3.0]]  
 [[1.0, 2.0], [4.0, 2.999]]

julia> isless0 = (x,y) -> [x[1][1],x[1][2],x[2][2]] < [y[1][1],y[1][2],y[2][2]]
#66 (generic function with 1 method)

julia> sort(lines, lt=isless0)
2-element Array{Array{Array{Float64,1},1},1}:
 [[1.0, 2.0], [4.0, 2.999]]
 [[1.0, 2.0], [3.0, 3.0]]  
```
"""
function presorted(V,EV, preview=false)
	lines = [[V[:,u],V[:,v]] for (u,v) in EV]
	orientedlines = Array{Array{Float64,1},1}[]
	for (v1,v2) in lines
		if v1 > v2 
			v1,v2 = v2,v1
		end
		push!(orientedlines, [v1,v2])
	end
	isless0 = (x,y) -> [x[1][1],x[1][2],x[2][2]] < [y[1][1],y[1][2],y[2][2]]
	inputlines = sort(orientedlines, lt=isless0)
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
		sizes = rand(Float64, 2)
		if grid==[1,1]
			V,(_,EV,_) = Lar.cuboid(corner,true,corner+sizes)
		else
			V,(_,EV,_) = Lar.cuboidGrid(grid, true)
			g = Lar.Struct([ Lar.t(corner...), Lar.s(sizes...), (V,EV) ])
			V,EV = Lar.struct2lar(g)
		end
		center = (corner + corner+sizes)/2
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
	theintersection(V::Lar.Points, EV::Lar.Cells)
		(e1::Int, e2::Int)::Union{Nothing, Array}
		
Compute the theintersection point of two line segments `line1` and `line2`,
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

julia> theintersection(line1,line2)
2-element Array{Float64,1}:
 0.08846469268137241
 0.9405292234774614 

julia> theintersection(line1,line2)
2-element Array{Float64,1}:
 0.08846469268137241
 0.9405292234774614 

julia> theintersection(line1,-1*line2) # returns a `Nothing` value

julia> 
```
"""
function theintersection(line1,line2)::Union{Nothing, Array}
	x1,y1,x2,y2 = vcat(line1...)
	x3,y3,x4,y4 = vcat(line2...)

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
			 # segments parallel: no theintersection
			 (β,α) = 999
		end
	end
	if 0<=α<=1 && 0<=β<=1
		p = [x1,y1] + α * ([x2,y2] - [x1,y1])
		return p
	end
end



"""
	swapsegments(SL,segA,segB)

Swap two segments at their theintersection event, by exchanging the
values of pairs (kay, value).
"""
function swapsegments(SL,segA,segB)
	# get segment keys in SL
	a = segA[1]; 
	b = segB[1]
	# semitokens of segA and segB in SL
	a_st = searchsortedfirst(SL,a)
	b_st = searchsortedfirst(SL,b)
	if a_st ≠ pastendsemitoken(SL) && b_st ≠ pastendsemitoken(SL)
		# recover values of old segments from SL
		k_a,v_a = deref((SL, a_st))
		k_b,v_b = deref((SL, b_st))
		# delete old segments from SL
		delete!((SL, a_st)) 
		delete!((SL, b_st))
		# swap the two segments
		# insert the new segments (sorted on keys) in SL
		st1 = insert!(SL, k_a, v_b)
		st2 = insert!(SL, k_b, v_a)
	end
	return st1,st2
end



"""
	selectsegmentneighbors(SL::SortedMultiDict,i)

Compute `segE1` and/or `segE2`, i.e. the one/two segments adjacent to `segE`, 
the segment of `E` event, on the `SL` sweepline. `segE2` should be *above* `segE`;
`segE2` should be *below* `segE`.
"""
function selectsegmentneighbor(SL,i)
	segE1 = []; segE2 = []
	if i ≠ pastendsemitoken(SL)
		st2 = advance((SL,i))
		keyE2 = st2 ≠ pastendsemitoken(SL) ? deref_key((SL,st2))[1] : nothing
		segE2 = keyE2 ≠ nothing ? deref_value((SL, st2)) : []
	end
	if i ≠ beforestartsemitoken(SL)
		st1 = regress((SL,i))
		keyE1 = st1 ≠ beforestartsemitoken(SL) ? deref_key((SL,st1))[1] : nothing
		segE1 = keyE1 ≠ nothing ? deref_value((SL, st1)) : []
	end
	return segE2, segE1
end


function sweepline(V,EV)
	# event generation and ordering
	segments = presortedlines(V,EV)
	evpairs = [[(v1,v2,"start",k), (v2,v1,"end",k)] for (k,(v1,v2)) in enumerate(segments)]
	events = sort(cat(evpairs))
	eventdict = Dict(zip([(e[4],e[3],e[1]) for e in events], events))

	# Initialize event queue ξ = all segment endpoints; Sort ξ by increasing x and y
	pqkeys = [(e[1],e[3],e[4]) for e in events]
	pqvalues = events
	ξ = PriorityQueue(zip(pqkeys,pqvalues))
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Any,Any}()
	# Initialized output theintersection list Λ to be empty
	Λ = Array{Float64,1}[]

	while length(ξ) ≠ 0 # (ξ is nonempty)

		obj = peek(ξ)
		E = obj[2]   # the  value of next event (k => v) from ξ 
		keyE = obj[1]   # the key of next event (k => v) from ξ 
		(v1, v2, nodetype, edgeid) = E
		# segE = E's segment
		e, segE = reverse(v1), (v1, v2, nodetype, edgeid) # key, value in SL
		if E[3] == "start" # (E is a left endpoint)
vals = [v for v in values(SL)]; for v in reverse(vals)	println(v) end
			# Add segE to SL
			i = insert!(SL, e, segE)  # SortedMultiDict, key, value -> semitoken
			if first(SL) !== last(SL) # more than one segment in SL	
				# compute segA and/or segB
				segA, segB = selectsegmentneighbor(SL,i)	
					
				if segA ≠ [] 
				# segA = the segment above segE in SL 
					a = segA[4]
					I = theintersection(segE,segA)
					# (if Intersect( segE with segA) exists)
					if typeof(I) ≠ Nothing
#						Insert I into ξ
#						key = (I,"int",edgeid); val = (I,I,"int",a)
#						enqueue!(ξ, key,val) 
					segB = segE
					# (I is not in ξ already) 
					# no problem in case I is overwritten
					akey = (segA[4],segA[3],segA[1])
					bkey = (segB[4],segB[3],segB[1])
					eventdict[akey] = (I,eventdict[segA][2],segA[3],segA[4])
					eventdict[bkey] = (I,eventdict[segB][2],segB[3],segB[4])
					# Insert I into ξ
					key = (I,"int",segA[4]); val = (I,I,"int",segB[4])
					enqueue!(ξ, key,val) 
					end
				end
				if segB ≠ [] 
				# segB = the segment below segE in SL 
					b = segB[4]
					I = theintersection(segE,segB)
					# (if Intersect( segE with segB) exists)
					if typeof(I) ≠ Nothing
#						Insert I into ξ
#						key = (I,"int",edgeid); val = (I,I,"int",b)
#						enqueue!(ξ, key,val) 
					segA = segE
					# (I is not in ξ already) 
					# no problem in case I is overwritten
					akey = (segA[4],segA[3],segA[1])
					bkey = (segB[4],segB[3],segB[1])
					eventdict[akey] = (I,eventdict[segA][2],segA[3],segA[4])
					eventdict[bkey] = (I,eventdict[segB][2],segB[3],segB[4])
					# Insert I into ξ
					key = (I,"int",segA[4]); val = (I,I,"int",segB[4])
					enqueue!(ξ, key,val) 
					end
				end
			end			
		elseif E[3] == "end" # (E is a right endpoint)
vals = [v for v in values(SL)]; for v in reverse(vals)	println(v) end
			# segE = E's segment
			# e = key(segE); st1 = corresponding semitoken in SL
			key = reverse(E[1])
			E_st = searchsortedlast(SL,key)
			if E_st ≠ beforestartsemitoken(SL) # no more events in SL 
				# compute segA and/or segB
				segA, segB = selectsegmentneighbor(SL,E_st)	
				#if segA==[] && segB==[] break end	
				# segA = the segment above segE in SL 
#				if segA ≠ []
#					a = segA[4] # edgeid of segA
#				end
#				# segB = the segment below segE in SL 
#				if segB ≠ []
#					b = segB[4]  # edgeid of segB
#				end
				# Remove segE from SL
				delete!((SL,E_st))
				# (I = Intersect( segA with segB) exists)
				if segA ≠ [] && segB ≠ []
					I = theintersection(segA, segB)
					if typeof(I) ≠ Nothing
						# (I is not in ξ already) 
						# no problem in case I is overwritten
						akey = (segA[4],segA[3],segA[1])
						bkey = (segB[4],segB[3],segB[1])
						eventdict[akey] = (I,eventdict[segA][2],segA[3],segA[4])
						eventdict[bkey] = (I,eventdict[segB][2],segB[3],segB[4])
						# Insert I into ξ
						key = (I,"int",segA[4]); val = (I,I,"int",segB[4])
						enqueue!(ξ, key,val) 
					end
				end
			end
		else # E is an theintersection event
#for (k,v) in SL @show (k,v) end
			# Add E to the output list Λ
			push!(Λ, E[1])
			# the two intersecting segments generating E
			a = eventdict[(keyE[3],"start")]
			b = eventdict[(E[4],"start")]
			a_st = searchsortedfirst(SL,reverse(a[1])) # KO: look via above/below
			b_st = searchsortedfirst(SL,reverse(b[1])) # KO: look via above/below
			segA = deref((SL,a_st))	
			segB = deref((SL,b_st))
			# Let segE1 above segE2 be E's intersecting segments in SL
			segE1,segE2 = compare(SL,a_st,b_st)==1 ? (segA,segB) : (segB,segA)
			# Swap their positions so that segE2 is now above segE1
			stE1,stE2 =  swapsegments(SL,segE1,segE2)
			# segA = the segment above segE2 in SL
			segA, segE2 = selectsegmentneighbor(SL,stE2)		
			# segB = the segment below segE1 in SL
			segE1, segB = selectsegmentneighbor(SL,stE1)					
			if segE1 ≠ [] && segB ≠ []
				I = theintersection(segE1,segB)
				# (I = Intersect(segE1 with segB) exists)
				if typeof(I) ≠ Nothing
					segA = segE1
					# (I is not in ξ already) 
					# no problem in case I is overwritten
					akey = (segA[4],segA[3],segA[1])
					bkey = (segB[4],segB[3],segB[1])
					eventdict[akey] = (I,eventdict[segA][2],segA[3],segA[4])
					eventdict[bkey] = (I,eventdict[segB][2],segB[3],segB[4])
					# Insert I into ξ
					key = (I,"int",segA[4]); val = (I,I,"int",segB[4])
					enqueue!(ξ, key,val) 
@show eventdict;
				end
			end
			if segA ≠ [] && segE2 ≠ []
				I = theintersection(segA, segE2)
				# (I = Intersect(segE2 with segA) exists)
				if typeof(I) ≠ Nothing
					segB = segE2
					# (I is not in ξ already) 
					# no problem in case I is overwritten
					akey = (segA[4],segA[3],segA[1])
					bkey = (segB[4],segB[3],segB[1])
					eventdict[akey] = (I,eventdict[segA][2],segA[3],segA[4])
					eventdict[bkey] = (I,eventdict[segB][2],segB[3],segB[4])
					# Insert I into ξ
					key = (I,"int",segA[4]); val = (I,I,"int",segB[4])
					enqueue!(ξ, key,val) 
				end
			end
			##stE1,stE2 = swapsegments(SL,segE1,segE2) ## ??? to remove
		end
		# remove E from ξ
		dequeue!(ξ)  
	end
	return Λ
end



# EXAMPLE 0
# data generation
lines = [[[1,3],[10,5]],[[2,6],[11,2.5]],
[[3,4],[8,2.25]],[[5,1.25],[12,5]]]
V,EV = lines2lar(lines)
Plasm.view(Plasm.numbering(2.)((V,[[[k] for k=1:size(V,2)], EV])))
# data sorting 
V = Plasm.normalize(V,flag=true)
W,EW = presorted(V,EV)
Plasm.view(Plasm.numbering(.25)((W,[[[k] for k=1:size(W,2)], EW ])))

Λ = sweepline(V,EV)


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

