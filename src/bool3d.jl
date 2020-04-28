using LinearAlgebraicRepresentation, SparseArrays
Lar = LinearAlgebraicRepresentation;
using IntervalTrees,LinearAlgebra
# using Revise,
# using OhMyREPL

#=
Method to compute an internal point to a polyhedron.
----------------------------------------------------

1. Take two of points close to the opposite sides of any face of a polyhedron, e.g., the first face.
2. For each of the two points compute the intersections of a (vertical) ray with the planes (of the faces) intersected by the ray (positive direction of the half-line).
3. Transform each such plane and face (and intersection point) to 2D space.
4. Test for point-in-polygon intersection.
5. Compute the parity of the intersection points for each ray.
6. Invariant:  if one is even; the other is odd.
7. The initial point with odd number of intersection points is interior to the polyhedron. The other is exterior.
=#

"""
	spaceindex(point3d)(model)

Compute the set of face boxes of possible intersection with a point-ray.
Work in 3D, where the ray direction is parallel to the z-axis.
Return an array of indices of face.

#	Example
```
julia> V,(VV,EV,FV,CV) = Lar.cuboidGrid([1,1,1],true)

julia> spaceindex([.5,.5,.5])((V,FV))
3-element Array{Int64,1}:
 5
 6
```
"""
function spaceindex(point3d::Array{Float64,1})::Function
	function spaceindex0(model::Lar.LAR)::Array{Int,1}
				V,CV = copy(model[1]),copy(model[2])
		V = [V point3d]
		dim, idx = size(V)
		push!(CV, [idx,idx,idx])
		cellpoints = [ V[:,CV[k]]::Lar.Points for k=1:length(CV) ]
		#----------------------------------------------------------
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		xboxdict = Lar.coordintervals(1,bboxes)
		yboxdict = Lar.coordintervals(2,bboxes)
		# xs,ys are IntervalTree type
		xs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in xboxdict
			xs[tuple(key...)] = boxset
		end
		ys = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in yboxdict
			ys[tuple(key...)] = boxset
		end
		xcovers = Lar.boxcovering(bboxes, 1, xs)
		ycovers = Lar.boxcovering(bboxes, 2, ys)
		covers = [intersect(pair...) for pair in zip(xcovers,ycovers)]

		# add new code part

		# remove each cell from its cover
		pointcover = setdiff(covers[end],[idx+1])
		return pointcover[1:end-1]
	end
	return spaceindex0
end

"""
	rayintersection(point3d::Array{Float64})(V,FV,face::Int)

Compute the intersection point of the vertical line through `point3d` w `face`.
If the face is parallel to `z axis` return `false`.
# Example
```
julia> V,(VV,EV,FV,CV) = Lar.simplex(3,true);

julia> V
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> FV
4-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 [1, 2, 4]
 [1, 3, 4]
 [2, 3, 4]

 julia> Lar.rayintersection([.333,.333,0])(V,FV,4)
 3-element Array{Float64,1}:
  0.333
  0.333
  0.3340000000000001
```
"""
function rayintersection(point3d)
	function rayintersection0(V, FV, face::Int)
						l0, l = point3d, [0,0,1.]
		ps = V[:,FV[face]]  # face points
		p0 = ps[:,1]
		v1, v2 = ps[:,2]-p0, ps[:,3]-p0
		n = LinearAlgebra.normalize(cross( v1,v2  ))

		denom = dot(n, l)
		if (abs(denom) > 1e-8) #1e-6
			p0l0 = p0 - l0
			t = dot(p0l0, n) / denom
			if t>0 return l0 + t*l end
		else
			#error("ray and face are parallel")
			return false
	 	end
	end
	return rayintersection0
end


"""
	planemap(V,copEV,copFE,face)(point)

Tranform the 3D face and the 3D point in their homologous 2D, in order to test for containment.
"""
function planemap(V,copEV,copFE,face)
	fv, edges = Lar.vcycle(copEV, copFE, face)
	# Fv = Dict(zip(1:length(fv),fv))
	# edges = [[Fv[u],Fv[v]] for (u,v) in edges]
	function planemap0(point)
		vs = V[:,fv]
		# Plasm.view(Plasm.numbering(0.5)((vs,[[[k] for k=1:4],edges])))
		#translation
		point = point .- vs[:,1]
		vs = vs .- vs[:,1]
		u,v = edges[1]
		z,w = [[z,w] for (z,w) in edges if z==v][1]
		v1 = vs[:,u]-vs[:,v]
		v2 = vs[:,w]-vs[:,v]
		v3 = cross(v2,v1)
		M = [v1 v2 v3]
		vs = inv(M) * [vs point]
		outvs = vs[1:2, 1:end-1]
		outpoint = vs[1:2, end]
		return outvs, edges, outpoint
	end
	return planemap0
end



function settestpoints(V,EV,FV,Fs, copEV,copFE)
	f = Fs[1]
	e = findnz(copFE[f,:])[1][1] # first (global) edge of first (global) face
	f1,f2 = findnz(copFE[:,e])[1] # two (global) faces incident on it
	v1,v2 = findnz(copEV[e,:])[1] # two (global) verts incident on in
	fdict = Dict(zip(Fs,1:length(Fs)))
	V1 = FV[fdict[f1]]
	V2 = FV[fdict[f2]]
	v1,v2 = intersect(V1,V2) # verified ... !
	t1 = V[:,v1], V[:,v2], V[:,[v for v in V1 if v≠v1 && v≠v2][1]]
	t2 = V[:,v2], V[:,v1], V[:,[v for v in V2 if v≠v1 && v≠v2][1]]
	n1 = LinearAlgebra.normalize(cross( t1[2]-t1[1], t1[3]-t1[1] ))
	n2 = LinearAlgebra.normalize(cross( t2[2]-t2[1], t2[3]-t2[1] ))
	p0 = (V[:,v1] + V[:,v2]) ./ 2
	n = n1 + n2
	ϵ = 1.0e-4
	ptest1 = p0 + ϵ * n
	ptest2 = p0 - ϵ * n
	return ptest1, ptest2
end


"""
	testinternalpoint(V::Lar.Points, EV::Lar.Cells, FV::Lar.Cells)

"""
function testinternalpoint(V,EV,FV)
	copEV = Lar.lar2cop(EV)
	copFV = Lar.lar2cop(FV)
	copFE = copFV * copEV'
	I,J,Val = findnz(copFE)
	triple = zip([(i,j,1) for (i,j,v) in zip(I,J,Val) if v==2]...)
	I,J,Val = map(collect,triple)
	Val = convert(Array{Int8,1},Val)
	copFE = sparse(I,J,Val)
	function testinternalpoint0(testpoint)
		intersectedfaces = Int64[]
		# spatial index for possible intersections with ray
		faces = spaceindex(testpoint)((V,FV))
		depot = []
		# face in faces :  indices of faces of possible intersection with ray
		for face in faces
			value = rayintersection(testpoint)(V,FV,face)
			if typeof(value) == Array{Float64,1} push!(depot, (face,value)) end
		end
		# actual containment test of ray point in faces within depot
		for (face,point3d) in depot
			vs, edges, point2d = planemap(V,copEV,copFE,face)(point3d)
			classify = Lar.pointInPolygonClassification(vs,edges)
			inOut = classify(point2d)
			if inOut!="p_out"
				push!(intersectedfaces,face)
			end
		end
		return intersectedfaces
	end
	return testinternalpoint0
end



"""
	getinternalpoint(V,EV,FV,Fs, copEV,copFE)

"""
function getinternalpoint(V,EV,FV,Fs, copEV,copFE)
	#edges for v1=FV[1][1]
	ptest1, ptest2 = Lar.settestpoints(V,EV,FV,Fs, copEV,copFE)
	intersectedfaces = Int64[]
	#GL.VIEW([ GL.GLFrame, GL.GLLines(V,EV), GL.GLPoints([ptest1'; ptest2']) ]);

	# for each test point compute the face planes intersected by vertical ray
	dep1, dep2 = [],[]
	# face in Fs : global indices of faces of current solid
	for (f,face) in enumerate(Fs)
		ret1 = rayintersection(ptest1)(V,FV,f)
		ret2 = rayintersection(ptest2)(V,FV,f)
		if typeof(ret1) == Array{Float64,1} push!(dep1, (face,ret1)) end
		if typeof(ret2) == Array{Float64,1} push!(dep2, (face,ret2)) end
	end

	# transform each plane in 2D and look whether the intersection point is internal
	# return the test point with odd numeber of ray intersections
	k1,k2 = 0,0
	for (face,point3d) in dep1
		vs, edges, point2d = planemap(V,copEV,copFE,face)(point3d)
		classify = Lar.pointInPolygonClassification(vs,edges)
		inOut = classify(point2d)
		if inOut!="p_out"
			k1+=1
			push!(intersectedfaces,face)
		end
	end
	if k1 % 2 == 1 return ptest1,intersectedfaces
	else
		for (face,point3d) in dep2
			vs, edges, point2d = planemap(V,copEV,copFE,face)(point3d)
			classify = Lar.pointInPolygonClassification(vs,edges)
			inOut = classify(point2d)
			if inOut!="p_out"
				k2+=1
				push!(intersectedfaces,face)
			end
		end
		if k2 % 2 == 1 return ptest2,intersectedfaces
		else println("error: while computing internal point of 3-cell") end
	end
end


# high level function

function chainbasis2solids(V,copEV,copFE,copCF)
	CF = [findnz(copCF[k,:])[1] for k=1:copCF.m]
	FE = [findnz(copFE[k,:])[1] for k=1:copFE.m]
	EV = [findnz(copEV[k,:])[1] for k=1:copEV.m]

	FEs = Array{Array{Int64,1},1}[]
	EVs = Array{Array{Array{Int64,1},1},1}[]
	FVs = Array{Array{Int64,1},1}[]
	for k=1:copCF.m
		push!( FEs, [collect(Set(cat([e for e in FE[f]]))) for f in CF[k]] )
		# edges in EVs are aggregated by face, in order to answer point-classifications
		push!( EVs, [[EV[e] for e in FE[f]] for f in CF[k]] )
		push!( FVs, [collect(Set(cat([EV[e] for e in FE[f]]))) for f in CF[k]] )
	end
	pols = collect(zip(EVs,FVs,FEs))
	W = convert(Lar.Points,V')
	return W,pols,CF
end



function internalpoints(V,copEV,copFE,copCF)
	# transform each 3-cell in a solid (via Lar model)
	#-------------------------------------------------------------------------------
	U,pols,CF = chainbasis2solids(V,copEV,copFE,copCF)
	# compute, for each `pol` (3-cell) in `pols`, one `internalpoint`.
	#-------------------------------------------------------------------------------
	innerpoints = []
	intersectedfaces = []
	for k=1:length(pols)
		(EV,FV,FE),Fs = pols[k],CF[k]
		EV = convert(Lar.Cells,collect(Set(cat(EV))))
		#GL.VIEW([ GL.GLFrame, GL.GLLines(V,EV) ]);
		points,facenumber = getinternalpoint(V,EV,FV,Fs, copEV,copFE)
		push!(innerpoints,points)
		push!(intersectedfaces,facenumber)
	end
	return innerpoints,intersectedfaces
end


################################################################################
#=
After the arrangement, extract all the d-cells from (d-1)-coboundary as isolated polyhedra.
Then compute a single interior point for each of them.
Then compare each such point against all input boundaries, in order to compute those which it was interior to. Extend this point membership as 3-cell containment within the relative input solids.
The point membership with a boundary consists in the parity count of the intersection points of a vertical ray starting at the test point, with the boundary surface.
=#
################################################################################

function bool3d(assembly)
	# input of affine assembly
	#-------------------------------------------------------------------------------
	V,FV,EV = Lar.struct2lar(assembly)
	cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
	cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells); ## TODO: debug
	W = convert(Lar.Points, V');
	# generate the 3D space arrangement
	#-------------------------------------------------------------------------------
	V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)
	W = convert(Lar.Points, V');
	#V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)
	innerpoints,intersectedfaces = Lar.internalpoints(W,copEV,copFE,copCF[2:end,:])
	# associate internal points to 3-cells
	#-------------------------------------------------------------------------------
	listOfModels = Lar.evalStruct(assembly)
	inputfacenumbers = [length(listOfModels[k][2]) for k=1:length(listOfModels)]
	cumulative = cumsum([0;inputfacenumbers]).+1
	fspans = collect(zip(cumulative[1:end-1], cumulative[2:end].-1))
	span(h) = [j for j=1:length(fspans) if fspans[j][1]<=h<=fspans[j][2] ]
	# test input data for containment of reference points
	#-------------------------------------------------------------------------------
	V,FV,EV = Lar.struct2lar(assembly)
	containmenttest = Lar.testinternalpoint(V,EV,FV)
	# currently copCF contains the outercell in first column ...
	# TODO remove first row and column, in case (look at file src/)
	boolmatrix = BitArray(undef, length(innerpoints)+1, length(fspans)+1)
	boolmatrix[1,1] = 1
	for (k,point) in enumerate(innerpoints) # k runs on columns
		cells = containmenttest(point) # contents of columns
		#println(k," ",faces)
		rows = [span(h) for h in cells]
		for l in cat(rows)
			boolmatrix[k+1,l+1] = 1
		end
	end
	return W, (copEV, copFE, copCF), boolmatrix
end
function bool3d(expr, assembly)
	V, copEV, copFE, copCF, boolmatrix = bool3d(assembly)
	# TODO please see template "bool2d()" in src/bool2d.jl
	return V, (copEV, copFE, copCF), boolmatrix, result
end

################################################################################

#
# #pointcover = spaceindex([0,0,0.])((V,FV))
