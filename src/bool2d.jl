using LinearAlgebraicRepresentation, SparseArrays
Lar = LinearAlgebraicRepresentation;
using IntervalTrees,LinearAlgebra
# using Revise, OhMyREPL


function settestpoints2d(W,EV,FV,f, copEV,copFE) # W by rows
	e = findnz(copFE[f,:])[1][1] # first edge of face f
	v1,v2 = findnz(copEV[e,:])[1] # two (global) verts incident on in
    t = W[v2,:] - W[v1,:]
    n = [-t[2],t[1]]
	p0 = (W[v1,:] + W[v2,:]) ./ 2
	ϵ = 1.0e-4
	ptest1 = p0 + ϵ * n
	ptest2 = p0 - ϵ * n
	return ptest1, ptest2
end


function getinternalpoint2d(W,EV,FV, f, copEV,copFE) # W by rows
	#edges for v1=FV[1][1]
	ptest1, ptest2 = Lar.settestpoints2d(W,EV,FV, f, copEV,copFE)
    edges = [findnz(copEV[e,:])[1] for e in findnz(copFE[f,:])[1]]
    V = convert(Lar.Points,W')
    classify = Lar.pointInPolygonClassification(V,edges)
    if classify(ptest1) == "p_in"
        return ptest1
    elseif classify(ptest2) == "p_in"
        return ptest2
    else
        error("classifying inner point in face $f")
    end
end



function chainbasis2polygons(V,copEV,copFE)
	FE = [findnz(copFE[k,:])[1] for k=1:copFE.m]
	EV = [findnz(copEV[k,:])[1] for k=1:copEV.m]

	FEs = Array{Int64,1}[]
	EVs = Array{Array{Int64,1},1}[]
	FVs = Array{Int64,1}[]
	for f=1:copFE.m
		push!( FEs, collect(Set(cat([e for e in FE[f]]))) )
		# edges in EVs are aggregated by face, in order to answer point-classifications
		push!( EVs, [EV[e] for e in FE[f]]  )
		push!( FVs, collect(Set(cat([EV[e] for e in FE[f]])))  )
	end
	polygons = collect(zip(EVs,FVs,FEs)) # (EV,FV,FFE)s for polygon
	W = convert(Lar.Points,V')
	return W,polygons,FE
end


function internalpoints2d(W,copEV,copFE) # W by rows
	# transform each 2-cell in a solid (via Lar model)
	#-------------------------------------------------------------------------------
	U,pols,FE = Lar.chainbasis2polygons(W,copEV,copFE) # V by rows
	# compute, for each `pol` (3-cell) in `pols`, one `internalpoint`.
	#-------------------------------------------------------------------------------
	internalpoints = []
	for f=1:length(pols)
		(EV,FV,FE) = pols[f]
		#GL.VIEW([ GL.GLFrame, GL.GLLines(V,EV) ]);
		internalpoint = Lar.getinternalpoint2d(W,EV,FV, f, copEV,copFE)
		push!(internalpoints,internalpoint)
	end
	return internalpoints
end

"""
	testinternalpoint2d(V::Lar.Points, EV::Lar.Cells, FV::Lar.Cells)

"""
function testinternalpoint2d(listOfModels)
	function testinternalpoint0(testpoint)
		intersectedfaces = Int64[]
		depot = []
		# actual containment test of ray point in faces within depot
		for (k,model) in enumerate(listOfModels)
            verts,edges = model
			classify = Lar.pointInPolygonClassification(verts,edges)
			inOut = classify(testpoint)
			if inOut == "p_in"
				push!(intersectedfaces,k)
			end
		end
		return intersectedfaces
	end
	return testinternalpoint0
end


################################################################################

function booleanops2d(assembly)
	# input of affine assembly
	#-------------------------------------------------------------------------------
	# V,EV = Lar.struct2lar(assembly) #TODO proper different method
	V,EV = Lar.struct2lar(assembly)
	cop_EW = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
	W = convert(Lar.Points, V');
	# generate the 3D space arrangement
	#-------------------------------------------------------------------------------
	W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
	#V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)
	innerpoints = Lar.internalpoints2d(W,copEV,copFE[1:end,:])
	# associate internal points to 2-cells
	#-------------------------------------------------------------------------------
	listOfModels = Lar.evalStruct(assembly)
	inputfacenumbers = [length(listOfModels[k][2]) for k=1:length(listOfModels)]
	# test input data for containment of reference points
	#-------------------------------------------------------------------------------
	booleanmatrix = zeros(Int8, length(innerpoints)+1, length(listOfModels)+1)
    booleanmatrix[1,1] = 1
    containmenttest = Lar.testinternalpoint2d(listOfModels)
	for (k,point) in enumerate(innerpoints) # k runs on columns
		cells = containmenttest(point) # contents of columns
		#println(k," ",faces)
		for l in cells
			booleanmatrix[k+1,l+1] = 1
		end
	end
	return booleanmatrix
end
