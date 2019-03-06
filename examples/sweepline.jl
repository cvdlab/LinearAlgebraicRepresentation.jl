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
			 # segments parallel: no intersection
			 (β,α) = 999
		end
	end
	if 0<=α<=1 && 0<=β<=1
		p = [x1,y1] + α * ([x2,y2] - [x1,y1])
		return p
	end
end


"""
	computesegposition(SL,h,xe,ye)

Compute the current position of a neighbourhood line with respect to a new segment.

`SL` (the sweep line) is a SortedMultiDict; `h` is the semitoken of a segment on SL;
`(xe, ye)` are the coordinates of vertex of a new segment lying on `SL`, w.r.t. 
compare the segment's `h` position.  The function returns a pair of values for `segA`
and `segB`, one of which is the empty array.
"""
function computesegposition(SL,h,xe,ye)
	segA = Array{Float64,1}[]; 
	segB = Array{Float64,1}[]
	u,v,type,edgeid = deref_value((SL,h))
	Dx = xe - u[1]
	Dy = Dx * (v[2]-u[2])/(v[1]-u[1])  ##### NO vertical edges allowed !!!
	y = u[2] + Dy
	if y > ye
		push!(segA, u, v, [edgeid])
	elseif y < ye
		push!(segB, u, v, [edgeid])
	elseif y == ye
		nothing::Nothing  # by now ... ;-)
	end
	return segA,segB
end


"""
	selectsegmentneighbors(SL::SortedMultiDict,i)

Compute `segA` and/or `segB`, i.e. the one/two segments adjacent to `segE`, 
the segment of `E` event, on the `SL` sweepline. `segA` should be *above* `segE`;
`segB` should be *below* `segE`.
"""
function selectsegmentneighbor(SL,i)
	seghA, seghB = Array{Float64,1}[], Array{Float64,1}[]
	segkA, segkB = Array{Float64,1}[], Array{Float64,1}[]

	h = regress((SL,i))
	if h ≠ beforestartsemitoken(SL)
		seghA,seghB = computesegposition(SL,h,xe,ye)
		k = regress((SL,h))
		if k ≠ beforestartsemitoken(SL)
			segkA,segkB = computesegposition(SL,k,xe,ye)
		end
	end
	segA = union(seghA, segkA); segB = union(seghB, segkB); 
	return segA, segB
end


function sweepline(V,EV)
	# event generation and ordering
	segments = presortedlines(V,EV)
	evpairs = [[(v1,v2,"start",k), (v2,v1,"end",k)] for (k,(v1,v2)) in enumerate(segments)]
	isless0 = (x,y) -> [x[1][1],x[1][2],x[2][2]] < [y[1][1],y[1][2],y[2][2]]
	events = sort(cat(evpairs),lt=isless0)
	# Initialize event queue ξ = all segment endpoints; Sort ξ by increasing x and y
	pqkeys = [(e[1],e[3],e[4]) for e in events]
	pqvalues = events
	ξ = PriorityQueue(zip(pqkeys,pqvalues))
@show ξ,0 println()
	# Proper ordering object for the `SortedMultiDic` used by sweep line SL
	##o = Lt((x,y) -> (x,y))	
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Any,Any}()##o)
@show SL, 0 println()
	# Initialized output intersection list Λ to be empty
	Λ = Array{Array{Float64,1},1}[]
@show Λ,0 println()

	while length(ξ) ≠ 0 # (ξ is nonempty)
		E = peek(ξ)[2]   # the next value of event (k => v) from ξ 
@show E, 1 println()
		(v1, v2, nodetype, edgeid) = E
		# segE = E's segment
		e, segE = v1, (v1, v2, nodetype, edgeid) # key, value
		if E[3] == "start" # (E is a left endpoint)
@show "start" println()
			# Add segE to SL
			i = insert!(SL, e, segE)  # SortedMultiDict, key, value -> semitoken
@show SL, 1 println()
@show i println()
			global xe,ye = e
			if first(SL) !== last(SL) # more than one segment in SL	
				# compute segA and/or segB
				segA, segB = selectsegmentneighbor(SL,i)		
				if segA ≠ [] 
				# segA = the segment above segE in SL 
					a = Int(segA[3][1]); segA = segA[1:2] # edgeid of segA
					I = intersection(segE,segA)
					# (if Intersect( segE with segA) exists)
					if typeof(I) ≠ Nothing
						# Insert I into ξ
						key = (I,"int",edgeid); val = (I,I,"int",a)
						enqueue!(ξ, key, val) 
@show ξ,6 println()
					end
				end
				if segB ≠ [] 
				# segB = the segment below segE in SL 
					b = Int(segB[3][1]); segB = segB[1:2] # edgeid of segB
					I = intersection(segE,segB)
					# (if Intersect( segE with segB) exists)
					if typeof(I) ≠ Nothing
						# Insert I into ξ
						key = (I,"int",edgeid); val = (I,I,"int",b)
						enqueue!(ξ, key, val) 
@show ξ,5 println()
					end
				end
			end			
		elseif E[3] == "end" # (E is a right endpoint)
@show "end" println()
			# segE = E's segment
			# e = key(segE); st1 = corresponding semitoken in SL
			st1,st2 = searchequalrange(SL::SortedMultiDict, e)
			# compute segA and/or segB
			segA, segB = selectsegmentneighbor(SL,st1)		
			# segA = the segment above segE in SL 
			if segA ≠ []
				a = Int(segA[3][1]) # edgeid of segA
				segA = segA[1:2] 
			end
			# segB = the segment below segE in SL 
			if segB ≠ []
				b = Int(segB[3][1])
				segB = segB[1:2] # edgeid of segB
			end
			# Remove segE from SL
@show st1,st2 println()
			#delete!((SL,st1))
println("############ ECCOMI !!!!!")
			# (I = Intersect( segA with segB) exists)
			if segA ≠ [] && segB ≠ []
				I = intersection(segA, segB)
				if typeof(I) ≠ Nothing
					# (I is not in ξ already) 
					# no problem anyway (in case I is overwritten)
					# Insert I into ξ
					key = (I,"int",a); val = (I,I,"int",b)
					enqueue!(ξ, key, val) 
@show ξ,4 println()
				end
			end
		else # E is an intersection event
		@assert E[3] == "int"
@show "int" println()
			# Add E to the output list Λ
			push!(Λ, [E[1],E[2]])
@show Λ,1 println()
			# segE1 above segE2 be E's intersecting segments in SL 
			i = insert!(SL, e, segE) ### ????
			segE1, segE2 = selectsegmentneighbor(SL,i)
			if segE1 ≠ []
				h = Int(segE1[3][1]) # index of segment segE1
				v1 = V[:,EV[h][1]] # key of segE1 in SL
				hst1,hst2 = searchequalrange(SL::SortedMultiDict, v1)
				# extreme semitokens of segE1 in SL
				# segB = the segment below segE1 in SL
@show SL,hst1 println()
				segE1, segB = selectsegmentneighbor(SL,hst1)		
				if segE1 ≠ [] && segB ≠ []
					I = intersection(segE1,segB)
					# (I = Intersect(segE1 with segB) exists)
					if typeof(I) ≠ Nothing
						# (if I is not in ξ already) no problem in case
						# Insert I into ξ
						key = (I,"int",edgeid); val = (I,I,"int",b) # edgeid ??
						enqueue!(ξ, key, val) 
@show ξ,5 println()
					end
				end
			end
			if segE2 ≠ []
				k = Int(segE2[3][1])  # index of segment segE2
				v1 = V[:,EV[k][1]]  # key of segE2 in SL
				kst1,kst2 = searchequalrange(SL::SortedMultiDict, v1)
				# extreme semitokens of segE2 in SL
				# Swap their positions so that segE2 is now above segE1  
				# ??? HOW TO DO IT ???
				# segA = the segment above segE2 in SL
@show SL,kst1 println()
				segA, segE2 = selectsegmentneighbor(SL,kst1)		
@show segA, segE2
				if segA ≠ [] && segE2 ≠ []
					I = intersection(segE2,segA)
					# (I = Intersect(segE2 with segA) exists)
					if typeof(I) ≠ Nothing
						# (if I is not in ξ already) no problem in case
						# Insert I into ξ
						key = (I,"int",edgeid); val = (I,I,"int",a) # edgeid ??
						enqueue!(ξ, key, val) 
@show ξ,4 println()
					end
				end
			end
			delete!((SL, i)) ### ????
		end
		# remove E from ξ
@show ξ,1 println()
		dequeue!(ξ)  
@show ξ,2 println()
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





#=
V = [0.435352 0.449902 0.658165 0.672715 0.76547 0.685644 0.542851 0.463025 0.438641 0.446501 0.406689 0.414549 0.145006 0.104172 0.218371 0.177538 0.129832 0.11199 0.151923 0.134081 0.206417 0.174285 0.208523 0.176391 0.988503 0.779401 0.904498 0.695396 0.0573599 0.232703 0.122291 0.297634 0.412384 0.42493 0.578685 0.59123 0.290235 0.0044308 0.285804 0.0; 0.290109 0.272658 0.475885 0.458435 0.274124 0.513271 0.199815 0.438962 0.813539 0.876549 0.817525 0.880535 0.947022 0.879789 0.902464 0.835231 0.416359 0.365524 0.408606 0.357771 0.345616 0.120288 0.345315 0.119988 0.286924 0.348144 0.0 0.0612204 0.416447 0.157502 0.460414 0.20147 0.759923 0.751232 1.0 0.99131 0.409458 0.428899 0.34432 0.363761]

EV = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [5, 7], [6, 8], [9, 10], [11, 12], [9, 11], [10, 12], [13, 14], [15, 16], [13, 15], [14, 16], [17, 18], [19, 20], [17, 19], [18, 20], [21, 22], [23, 24], [21, 23], [22, 24], [25, 26], [27, 28], [25, 27], [26, 28], [29, 30], [31, 32], [29, 31], [30, 32], [33, 34], [35, 36], [33, 35], [34, 36], [37, 38], [39, 40], [37, 39], [38, 40]]

sweepline(V,EV)
=#