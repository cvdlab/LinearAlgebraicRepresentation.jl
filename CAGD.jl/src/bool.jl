function bool3(model::CAGD.Model)
    arranged_model = CAGD.spatial_arrangement(model, outerCell = true)
    innerpoints, intersectedfaces = CAGD.internalpoints(arranged_model)
    # number of input complexes
    nInC = size(model, 3, 1)
    infaces = [CAGD.getModelLoCell(model, 3, i) for i = 1 : nInC]
    # Determine in which initial element was h
    span(h) = [j for j = 1 : nInC if h in infaces[j]]
	containmenttest = CAGD.testinternalpoint(model)
	# currently copCF contains the outercell in first column ...
	# TODO remove first row and column, in case (look at file src/)
	boolmatrix = BitArray(undef, length(innerpoints)+1, nInC+1)
	boolmatrix[1,1] = 1
	for (k, point) in enumerate(innerpoints) # k runs on columns
		cells = containmenttest(point) # contents of columns
		#println(k," ",faces)
		rows = [span(h) for h in cells]
		for l in cat(rows)
			boolmatrix[k+1, l+1] = 1
		end
	end
	return arranged_model, boolmatrix
end

function internalpoints(model::CAGD.Model)
    return Lar.internalpoints(model.G, model.T[1], model.T[2], model.T[3][2:end, :])
end

"""
	testinternalpoint(V::Lar.Points, EV::Lar.Cells, FV::Lar.Cells)

"""
function testinternalpoint(model)
	copEV = abs.(model.T[1])
    copFE = abs.(model.T[2])
    V = model.G
    EV = Lar.cop2lar(model.T[1])
    FV = Lar.cop2lar(map(x -> floor(Int8, x / 2), copFE * copEV))
    
	I,J,Val = findnz(copFE)
	Val = convert(Array{Int8,1},Val)
	copFE = sparse(I,J,Val)
	function testinternalpoint0(testpoint)
		intersectedfaces = Int64[]
		# spatial index for possible intersections with ray
		faces = Lar.spaceindex(testpoint)((V, FV))
		depot = []
		# face in faces :  indices of faces of possible intersection with ray
		for face in faces
			value = Lar.rayintersection(testpoint)(V, FV, face)
			if typeof(value) == Array{Float64,1} push!(depot, (face,value)) end
		end
		# actual containment test of ray point in faces within depot
		for (face,point3d) in depot
			vs, edges, point2d = Lar.planemap(V, copEV, copFE, face)(point3d)
			classify = Lar.pointInPolygonClassification(vs,edges)
			inOut = classify(point2d)
			if inOut != "p_out"
				push!(intersectedfaces,face)
			end
		end
		return intersectedfaces
	end
	return testinternalpoint0
end
