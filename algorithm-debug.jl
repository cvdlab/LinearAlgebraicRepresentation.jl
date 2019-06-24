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
	# Proper ordering object for the `SortedMultiDic` used by sweep line SL
	##o = Lt((x,y) -> (x,y))
	# Initialize sweep line SL to be empty
	SL = SortedMultiDict{Any,Any}()##o)
	# Initialized output intersection list Λ to be empty
	Λ = Array{Array{Float64,1},1}[]

	while length(ξ) ≠ 0 # (ξ is nonempty)
		E = peek(ξ)[2]   # the next value of event (k => v) from ξ
		(v1, v2, nodetype, edgeid) = E
		# segE = E's segment
		e, segE = v1, (v1, v2, nodetype, edgeid) # key, value
		if E[3] == "start" # (E is a left endpoint)
			# Add segE to SL
			i = insert!(SL, e, segE)  # SortedMultiDict, key, value -> semitoken
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
						#(keyE,newsegE),(keyA,newsegA) = swapsegments(SL,I,segE,segA)
						key = (I,"int",edgeid); val = (I,I,"int",a)
						enqueue!(ξ, key,val)
						#enqueue!(ξ, keyA,newsegA)
					end
				end
				if segB ≠ []
				# segB = the segment below segE in SL
					b = Int(segB[3][1]); segB = segB[1:2] # edgeid of segB
					I = intersection(segE,segB)
					# (if Intersect( segE with segB) exists)
					if typeof(I) ≠ Nothing
						# Insert I into ξ
						#(keyE,newsegE),(keyB,newsegB) = swapsegments(SL,I,segE,segB)
						key = (I,"int",edgeid); val = (I,I,"int",b)
						enqueue!(ξ, key,val)
						#enqueue!(ξ, keyB,newsegB)
					end
				end
			end
		elseif E[3] == "end" # (E is a right endpoint)
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
			##delete!((SL,st1))
			# (I = Intersect( segA with segB) exists)
			if segA ≠ [] && segB ≠ []
				I = intersection(segA, segB)
				if typeof(I) ≠ Nothing
					# (I is not in ξ already)
					# no problem anyway (in case I is overwritten)
					# Insert I into ξ
					#(keyA,newsegA),(keyB,newsegB) = swapsegments(SL,I,segA,segB)
					key = (I,"int",a); val = (I,I,"int",b)
					enqueue!(ξ, key,val)
					#enqueue!(ξ, keyB,newsegB)
				end
			end
		else # E is an intersection event
		@assert E[3] == "int"
			# Add E to the output list Λ
			push!(Λ, [E[1],E[2]])
			# segE1 above segE2 be E's intersecting segments in SL
			i = insert!(SL, e, segE) ### ????
			segE1, segE2 = selectsegmentneighbor(SL,i)
			if segE1 ≠ []
				h = Int(segE1[3][1]) # index of segment segE1
				v1 = V[:,EV[h][1]] # key of segE1 in SL
				hst1,hst2 = searchequalrange(SL::SortedMultiDict, v1)
				# extreme semitokens of segE1 in SL
				# segB = the segment below segE1 in SL
				segE1, segB = selectsegmentneighbor(SL,hst1)
				if segE1 ≠ [] && segB ≠ []
					I = intersection(segE1,segB)
					# (I = Intersect(segE1 with segB) exists)
					if typeof(I) ≠ Nothing
						# (if I is not in ξ already) no problem in case
						# Insert I into ξ
						(keyE1,newsegE1),(keyB,newsegB) = swapsegments(SL,I,segE1,segB)
						key = (I,"int",edgeid); val = (I,I,"int",b) # edgeid ??
						enqueue!(ξ, key,val)
						#enqueue!(ξ, keyB,newsegB)
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
				segA, segE2 = selectsegmentneighbor(SL,kst1)
				if segA ≠ [] && segE2 ≠ []
					I = intersection(segA, segE2)
					# (I = Intersect(segE2 with segA) exists)
					if typeof(I) ≠ Nothing
						# (if I is not in ξ already) no problem in case
						# Insert I into ξ
						(keyA,newsegA),(keyE2,newsegE2) = swapsegments(SL,I,segA,segE2)
						##key = (I,"int",edgeid); val = (I,I,"int",a) # edgeid ??
						enqueue!(ξ, keyA,newsegA)
						enqueue!(ξ, keyE2,newsegE2)
					end
				end
			end
			delete!((SL, i)) ### ????
		end
		# remove E from ξ
		dequeue!(ξ)
	end
	return Λ
end
