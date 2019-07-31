using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using IntervalTrees,LinearAlgebra
using Revise, OhMyREPL

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
#@show point3d
	function spaceindex0(model::Lar.LAR)::Array{Int,1}
		#@show model
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

 julia> rayintersection([.333,.333,0])(V,FV,4)
 3-element Array{Float64,1}:
  0.333
  0.333
  0.3340000000000001
```
"""
function rayintersection(point3d)
#@show point3d
	function rayintersection0(V, FV, face::Int)
		#@show V;
	#@show FV;
	#@show face
		l0, l = point3d, [0,0,1.]
		ps = V[:,FV[face]]  # face points
		p0 = ps[:,1]
		v1, v2 = ps[:,2]-p0, ps[:,3]-p0
		n = normalize(cross( v1,v2  ))

		denom = dot(n, l)
		if (abs(denom) > 1e-6)
			p0l0 = p0 - l0
			t = dot(p0l0, n) / denom
			if t>0 return l0 + t*l end
		else
			#error("ray and face are parallel")
			return false
	 	end
		return rayintersection0
	end
end


"""
	removeconstrow(A::Array{Float64,2})::

Remove a row of constant values from a matrix.
"""
function removeconstrow(A::Array{Float64,2})
#@show A
	B = Array{Float64,1}[]
	global h = 0
	for k=1:size(A,1)
		rowtest = [el for el in A[k,:] if abs(el-A[k,1])<1e-3]
		if length(rowtest)!=0
			push!(B,A[k,:])
		else
			h = k
		end
	end
	B = convert(Array{Float64,2},hcat(B...)')
	return B,h
end


"""
	planemap(V,copEV,copFE,face)(point)

Tranform the 3D face and the 3D point in their homologous 2D, in order to test for containment.
"""
function planemap(V,copEV,copFE,face)
#@show findnz(copEV);
#@show findnz(copFE);
#@show face
	fv, edges = Lar.vcycle(copEV, copFE, face)
@show fv
@show edges
	function planemap0(point)
		#@show point
		vs = V[:,fv]
@show vs
		vs,h = removeconstrow(vs)
@show vs,h
		if h==0
			u,v = edges[1]
			z,w = [[z,w] for (z,w) in edges if z==v][1]
			v1 = vs[:,u]-vs[:,v]
			v2 = vs[:,w]-vs[:,v]
			v3 = cross(v2,v1)
			M = [v1 v2 v3]
@show M
			vs = inv(M) * vs
			vs,h = removeconstrow(vs)
		end

		return vs, edges, [point[k] for k=1:length(point) if k≠h] # TODO: debug
	end
	return planemap0
end


"""
	getinternalpoint(V::Lar.Points, FV::Lar.Cells)::Array(Float64)

"""
function getinternalpoint(V,EV,FV,Fs, copEV,copFE)
#@show V;
#@show FV;
#@show findnz(copEV);
#@show findnz(copFE);
	# get two test points close to the two sides of first face
	(v1,v2),v3 = EV[1],[v for (u,v) in EV if u==EV[1][2]][1]
	ps = [V[:,v1] V[:,v2] V[:,v3]]  # face points
	p0 = (ps[:,1]+ps[:,2])./2
	#GL.VIEW([ GL.GLFrame, GL.GLLines(U,EV), GL.GLPoints([ps p0]) ]);
	t = ps[:,2]-ps[:,1] # suppose first 3 points not aligned
	v1,v2 = ps[:,2]-ps[:,1], ps[:,3]-ps[:,2]
	n = normalize(cross( t, v2  ))
	ϵ = 1.0e-2
	ptest1 = p0 + ϵ*v1 + ϵ*v2 + ϵ*n  # point test one
	ptest2 = p0 + ϵ*v1 + ϵ*v2 - ϵ*n  # point test two
	GL.VIEW([ GL.GLFrame, GL.GLLines(U,EV), GL.GLPoints([ptest1'; ptest2']) ]);

	# for each test point compute the face planes intersected by vertical ray
	dep1, dep2 = [],[]
	# face in Fs : global indices of faces of current solid
	for (f,face) in enumerate(Fs)
		ret1 = rayintersection(ptest1)(V,FV,f)
		ret2 = rayintersection(ptest2)(V,FV,f)
		if typeof(ret1) == Array{Float64,1} push!(dep1, (face,ret1)) end
		if typeof(ret2) == Array{Float64,1} push!(dep2, (face,ret2)) end
	end
	p_on_planes = hcat([ret1 for (face,ret1) in dep1]...)
	fv, edges = Lar.vcycle(copEV, copFE, Fs[1])
	GL.VIEW([GL.GLFrame, GL.GLLines(V,EV),
		GL.GLPoints(convert(Array{Float64,2},p_on_planes')) ]);


	# transform each plane in 2D and look whether the intersection point is internal
	# return the test point with odd numeber of ray intersections
	k1,k2 = 0,0
	for (face,point3d) in dep1
		vs, edges, point2d = planemap(V,copEV,copFE,face)(point3d)
# p = convert(Array{Float64,2}, point2d')
# GL.VIEW([GL.GLFrame2, GL.GLLines(vs,edges), GL.GLPoints(p)])
		classify = Lar.pointInPolygonClassification(vs,edges)
		inOut = classify(point2d)
		println(inOut)
		if inOut!="p_out"  k1+=1 end
	end
	if k1 % 2 == 1 return ptest1
	else
		for (face,point3d) in dep2
			vs, edges, point2d = planemap(V,copEV,copFE,face)(point3d)
			classify = Lar.pointInPolygonClassification(vs,edges)
			inOut = classify(point2d)
			println(inOut)
			if inOut!="p_out"  k2+=1 end
		end
		if k2 % 2 == 1 return ptest2
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

################################################################################
#=
After the arrangement, extract all the d-cells from (d-1)-coboundary as isolated polyhedra.
Then compute a single interior point for each of them.

Then compare each such point against all input boundaries, in order to compute those which it was interior to. Extend this point membership as 3-cell containment within the relative input solids.

The point membership with a boundary consists in the parity count of the intersection points of a vertical ray starting at the test point, with the boundary surface.
=#
################################################################################

# Example generation
#-------------------------------------------------------------------------------
n,m,p = 1,1,1
V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
cube = V,FV,EV

threecubes = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube ])
V,FV,EV = Lar.struct2lar(threecubes)
GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame ]);
cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

# Arrangement computation
#-------------------------------------------------------------------------------
# generate the 3D space arrangement
V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)

# transform each 3-cell in a solid (via Lar model)
#-------------------------------------------------------------------------------
U,pols,CF = chainbasis2solids(V,copEV,copFE,copCF)

# compute, for each `pol` (3-cell) in `pols`, one `internalpoint`.
#-------------------------------------------------------------------------------


internalpoints = []
for k=1:length(pols)
@show k
	(EV,FV,FE),Fs = pols[k],CF[k]
	EV = convert(Lar.Cells,cat(EV))
	#GL.VIEW([ GL.GLFrame, GL.GLLines(U,EV) ]);
	internalpoint = getinternalpoint(U,EV,FV,Fs, copEV,copFE)
	push!(internalpoints,internalpoint)
end


#
#  points = [[-0.001, 0.898294, 0.14675],
#  [0.001, 0.898294, 0.14675] ,
#  [0.001, 0.898294, 0.14675] ,
#  [-0.001, 0.898294, 0.14675]]
#
# points = convert(Array{Float64,2}, hcat(points...)')
# V,FV,EV = Lar.struct2lar(threecubes)
# GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame, GL.GLPoints(points) ]);
#
#
# #pointcover = spaceindex([0,0,0.])((V,FV))
