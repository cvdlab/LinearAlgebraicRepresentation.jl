using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using IntervalTrees,LinearAlgebra

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
function spaceindex(point3d::Array{Float64,1},testpoint=true)::Function
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

 julia> rayintersection([.333,.333,0])(V,FV,4)
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
		n = normalize(cross( v1,v2  ))

		denom = dot(n, l)
		if (abs(denom) > 1e-6)
			p0l0 = p0 - l0
			t = dot(p0l0, n) / denom
			return l0 + t*l
		else
			#error("ray and face are parallel")
			return false
	 	end
		return rayintersection0
	end
end


function removeconstrow(A::Array{Float64,2})
	B = Array{Float64,1}[]
	global h = 0
	for k=1:size(A,1)
		rowtest = [el for el in A[k,:] if el≠A[k,1]]
		if length(rowtest)!=0
			push!(B,A[k,:])
		else
			h = k
		end
	end
	return B,h
end


"""
	planemap(V,copEV,copFE,face)(point)

Tranform the 3D face and the 3D point in their homologous 2D, in order to test for containment.
"""
function planemap(V,copEV,copFE,face)
	fv, edges = Lar.vcycle(copEV, copFE, face)
	function planemap0(point)
		vs = V[:,fv]
		vs,h = removeconstrow(vs)
		if h==0
			u,v = edges[1]
			z,w = [[z,w] for (z,w) in edges if z==v][1]
			v1 = V[:,u]-V[:,v]
			v2 = V[:,w]-V[:,v]
			v3 = cross(v2,v1)
			M = [v1 v2 v3]
			vs = (inv(M)*vs)
			vs,h = removeconstrow(vs)
		end
		return vs, [point[k] for k=1:length(point) if k≠h]
	end
	return planemap0
end


"""
	getinternalpoint(V::Lar.Points, FV::Lar.Cells)::Array(Float64)

"""
function getinternalpoint(V,FV)
	# get two test points close to the two sides of any face (first is OK)
	ps = V[:,FV[1]]  # face points
	p0 = ps[:,1]
	v1, v2 = ps[:,2]-p0, ps[:,3]-p0 # suppose first 3 points not aligned
	n = normalize(cross( v1,v2  ))
	ϵ = 1.0e-3
	ptest1 = p0 + ϵ*v1 + ϵ*v2 + ϵ*n  # point test one
	ptest2 = p0 + ϵ*v1 + ϵ*v2 - ϵ*n  # point test two

	# for each test point compute the face planes intersected by vertical ray
	dep1, dep2 = [],[]
	for face in 1:length(FV)
		ret1 = rayintersection(ptest1)(V,FV,face)
		ret2 = rayintersection(ptest2)(V,FV,face)
		if ret1 ≠ false push!(dep1, ret1) end
		if ret2 ≠ false push!(dep2, ret2) end
	end

	# transform each plane in 2D and look whether the intersection point is internal


	# return the test point with even numeber of intersections

end

# Example generation


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

pointcover = Lar.spaceindex([0,0,0.])((V,FV))

VV = [[k] for k=1:size(V,2)]
GL.VIEW( push!(GL.numbering(0.5)( (V,[ VV,EV,FV ]), GL.COLORS[1],1 ), GL.GLFrame) )



V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)



W = convert(Lar.Points, V')
