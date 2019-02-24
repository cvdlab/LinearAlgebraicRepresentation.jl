using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
using IntervalTrees
using SparseArrays



#---------------------------------------------------------------------
#	2D containment test
#---------------------------------------------------------------------

""" 

Half-line crossing test. Utility function for `pointInPolygonClassification` function.
Update the `count` depending of the actual crossing of the tile half-line.
"""
function crossingTest(new::Int, old::Int, count::Float64, status::Int)::Number
    if status == 0
        status = new
        count += 0.5
    else
        if status == old
        	count += 0.5
        else
        	count -= 0.5
        end
        status = 0
    end
end



""" 
	setTile(box)(point)

Set the `tileCode` of the 2D bbox `[b1,b2,b3,b4]:=[ymax,ymin,xmax,xmin]:= x,x,y,y` 
including the 2D `point` of `x,y` coordinates.
Depending on `point` position, `tileCode` ranges in `0:15`, and uses bit operators.
Used to set the plane tiling depending on position of the query point, 
in order to subsequently test the tile codes of edges of a 2D polygon, and determine 
if the query point is either internal, external, or on the boundary of the polygon.
Function to be parallelized ...

```julia
c1,c2 = tilecode(p1),tilecode(p2)
c_edge, c_un, c_int = c1 ⊻ c2, c1 | c2, c1 & c2
```

"""
function setTile(box)
	tiles = [[9,1,5],[8,0,4],[10,2,6]]
	b1,b2,b3,b4 = box
	function tileCode(point)
		x,y = point
		code = 0
		if y>b1 code=code|1 end
		if y<b2 code=code|2 end
		if x>b3 code=code|4 end
		if x<b4 code=code|8 end
		return code
	end
	return tileCode
end



""" 
	pointInPolygonClassification(V,EV)(pnt)

Point in polygon classification.

# Example

```julia
result = []
classify = pointInPolygonClassification(V,EV)
```
"""
function pointInPolygonClassification(V,EV)
    function pointInPolygonClassification0(pnt)
        x,y = pnt
        xmin,xmax,ymin,ymax = x,x,y,y
        tilecode = setTile([ymax,ymin,xmax,xmin])
        count,status = 0,0
    
        for (k,edge) in enumerate(EV)
            p1,p2 = V[:,edge[1]],V[:,edge[2]]
            (x1,y1),(x2,y2) = p1,p2
            c1,c2 = tilecode(p1),tilecode(p2)
            c_edge, c_un, c_int = c1⊻c2, c1|c2, c1&c2
            
            if (c_edge == 0) & (c_un == 0) return "p_on" 
            elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
            elseif c_edge == 3
                if c_int == 0 return "p_on"
                elseif c_int == 4 count += 1 end
            elseif c_edge == 15
                x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
                if x_int > x count += 1
                elseif x_int == x return "p_on" end
            elseif (c_edge == 13) & ((c1==4) | (c2==4))
                    crossingTest(1,2,status,count)
            elseif (c_edge == 14) & ((c1==4) | (c2==4))
                    crossingTest(2,1,status,count)
            elseif c_edge == 7 count += 1
            elseif c_edge == 11 count = count
            elseif c_edge == 1
                if c_int == 0 return "p_on"
                elseif c_int == 4 crossingTest(1,2,status,count) end
            elseif c_edge == 2
                if c_int == 0 return "p_on"
                elseif c_int == 4 crossingTest(2,1,status,count) end
            elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
            elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
            elseif c_edge == 5
                if (c1==0) | (c2==0) return "p_on"
                else crossingTest(1,2,status,count) end
            elseif c_edge == 6
                if (c1==0) | (c2==0) return "p_on"
                else crossingTest(2,1,status,count) end
            elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
            elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
            end
        end
        if (round(count)%2)==1 
        	return "p_in"
        else 
        	return "p_out"
        end
    end
    return pointInPolygonClassification0
end



#---------------------------------------------------------------------
#	Refactoring pipeline
#---------------------------------------------------------------------

"""
	input_collection(data::Array)::Tuple

*Facet selection*. Construction of a ``(d-1)``-dimensional collection from a ``(d-1)``- 
or ``d``-dimensional one. ``0-chain`` of `LAR` type are used as *input*.

*Output* is ``admissible input`` for algorithms of the *2D/3D arrangement* pipeline.

# Example 2D

An assembly of geometric objects is generated, and their assembly, including rotated 
and translated chains, is built producing a collection of input LAR models.

```julia
V,(_,EV,FV) = Lar.cuboidGrid([4,4],true);
W,(_,EW,FW) = Lar.cuboidGrid([3,5],true);
mycircle(r,n) = Lar.circle(r)(n)

data2d1 = (V,EV)
data2d2 = Lar.Struct([ Lar.t(2,2), Lar.r(pi/3), Lar.t(-1.5,-2.5), (W,EW) ])
data2d3 = Lar.Struct([ Lar.t(2,2), mycircle(2.5,16) ])
data2d4 = Lar.Struct([ Lar.t(3.5,3.5), mycircle(.25,16) ])
data2d5 = Lar.Struct([ Lar.t(5,3.5), mycircle(.5,16) ])
data2d6 = Lar.Struct([ Lar.t(5,3.5), mycircle(.25,16) ])

model2d = input_collection( [ data2d1, data2d2, data2d3, data2d4, data2d5, data2d6 ] )
V,EV = model2d
VV = [[k] for k in 1:size(V,2)];
Plasm.view( Plasm.numbering(.5)((V,[VV,EV])) )
```
Note that `V,EV` is not a cellular complex, since 1-cells intersect out of 0-cells.

# Example 3D

```julia
V,FV = Lar.sphere(2)([3,4])
EV = Lar.simplexFacets(FV)
mysphere = V,FV,EV

data3d1 = mysphere
data3d2 = Lar.Struct([ Lar.t(0,1,0), mysphere ])
data3d3 = Lar.Struct([ Lar.t(0,0.5,0), Lar.s(0.4,0.4,0.4), mysphere ])
data3d4 = Lar.Struct([ Lar.t(4,0,0), Lar.s(0.8,0.8,0.8), mysphere ])
data3d5 = Lar.Struct([ Lar.t(4,0,0), Lar.s(0.4,0.4,0.4), mysphere ])

model3d = input_collection([ data3d1, data3d2, data3d3, data3d4, data3d5 ])
V,FV,EV = model3d
VV = [[k] for k in 1:size(V,2)];
Plasm.view( Plasm.numbering(.6)((V,[VV, EV])) )
```

Note that `V,FV,EV` is not a cellular complex, since 1-cells and
2-cells intersect out of 0-cells.

"""
function input_collection(data::Array)::Lar.LAR
	assembly = Lar.Struct(data)
	return Lar.struct2lar(assembly)
end




"""
	indexing(model::Lar.LAR)
	
*Spatial index* made by ``d`` *interval-trees* on 
bounding boxes of ``sigma in S_{d−1}``. Spatial queries solved by
intersection of ``d`` queries on IntervalTrees generated by
bounding-boxes of geometric objects (LAR cells).

# Example 2D

```julia
model = model2d
xs,ys =  indexing(model2d);
```

# Example 3D

```julia
model = model3d
xs,ys,zs =  indexing(model3d);
```
"""
function indexing(model)
	V,CV = model[1:2]
	dim = size(V,1)
	@assert length(model) == dim  #n. chains == dim space
	
	cellpoints = [ V[:,CV[k]]::Lar.Points for k=1:length(CV) ]
	bboxes = [Lar.bbox(cell) for cell in cellpoints]
	ordboxes = [map(x->sort(x,dims=2),bboxes[k]) for k=1:length(bboxes)]
	xordboxes = sort([(vec(ordboxes[h][1]), h) for h=1:length(ordboxes)])
	yordboxes = sort([(vec(ordboxes[h][2]), h) for h=1:length(ordboxes)])
	
	# xs IntervalTree
	xintervals = IntervalValue{Float64,Int}[]
	for (k,xbox) in enumerate(xordboxes)
		push!(xintervals, IntervalValue{Float64,Int}(xbox[1][1], xbox[1][2], xbox[2]))
	end
	xs = IntervalTree{Float64, IntervalValue{Float64,Int} }(xintervals)
	
	# ys IntervalTree
	yintervals = IntervalValue{Float64,Int}[]
	for (k,ybox) in enumerate(yordboxes)
		push!(yintervals, IntervalValue{Float64,Int}(ybox[1][1], ybox[1][2], ybox[2]))
	end
	ys = IntervalTree{Float64, IntervalValue{Float64,Int}}(yintervals)
	
	if dim == 3
		# zs IntervalTree
		zintervals = IntervalValue{Float64,Int}[]
		for (k,zbox) in enumerate(zordboxes)
			push!(zintervals, IntervalValue{Float64,Int}(zbox[1][1], zbox[1][2], zbox[2]))
		end
		zs = IntervalTree{Float64, IntervalValue{Float64,Int}}(zintervals)
	end
	if dim==2 return (xs,ys)
	elseif dim==3 return (xs,ys,zs) 
	end
end



"""
	spaceindex(model::Lar.LAR)::Array{Array{Int,1},1}
	
Generation of *space indexes* for all ``(d-1)``-dim cell members of `model`.
`model` input must be a pair for 2d LAR, and a triple for 3d LAR.
The return value is an array of arrays of `int`s, indexing 2-, 3-cells whose 
containment boxes are intersecting the containment box of the first cell. 
According to Hoffmann, Hopcroft, and Karasick (1989) the worst-case complexity of
Boolean ops on such complexes equates the total sum of such numbers. 

# Example 2D

```julia
model = model2d
Sigma =  spaceindex(model2d);
typeof(Sigma) = Array{Array{Int,1},1}
```

# Example 3D

```julia
model = model3d
Sigma =  spaceindex(model3d);
typeof(Sigma) = Array{Array{Int,1},1}
```
"""
function spaceindex(model::Lar.LAR)::Array{Array{Int,1},1}
	V,CV = model[1:2]
	dim = size(V,1)
	@assert length(model) == dim  #n. chains == dim space
	
	if dim == 3 
		xs,ys,zs = Lar.indexing(model)
	elseif dim == 2 
		xs,ys = Lar.indexing(model)
	end
	
	function bbox(vertices::Lar.Points)
	   minimum = mapslices(x->min(x...), vertices, dims=2)
	   maximum = mapslices(x->max(x...), vertices, dims=2)
	   return minimum, maximum
	end
	
	spatialindex = []
	for (k,sigma) in enumerate(CV)
		Sigma = [k]
		cellpoints = V[:,sigma]
		vmin, vmax = map(x->sort(x,dims=2),bbox(cellpoints)) 		
		xquery = intersect( xs::IntervalTree, (vmin[1], vmax[1]) )
		xqs = [xint.value for xint in xquery]
		yquery = intersect( ys::IntervalTree, (vmin[2], vmax[2]) )
		yqs = [yint.value for yint in yquery]
		
		if dim == 3 
			zquery = intersect( zs::IntervalTree, (vmin[3], vmax[3]) )
			zqs = [zint.value for zint in zquery]
			xyzs = intersect(xqs,yqs,zqs)
		elseif dim == 2 
			xyzs = intersect(xqs,yqs)
		end
		
		I_sigma = xyzs #filter(x -> x>k,sort(xyzs))
		append!(Sigma, I_sigma)		
		push!(spatialindex,Sigma)
	end
	return spatialindex
end



"""
	decomposition2d(model::Lar.LAR)
	
Test function Lar.Arrangement.planar_arrangement.
Pairwise *intersection* of 2D *line segments* in ``σ ∪ I(σ)``, for each ``σ ∈ Sd−1``.

# Example 2D

```julia
V,EV = model2d
W, copEW, copFE = Lar.decomposition2d(model2d) # OK
Plasm.viewexploded(W, copEW)(1.2,1.2,1.2)
```
"""
function decomposition2d(model::Lar.LAR)
	V,EV = model
	spatialindex = spaceindex(model)
	
#test previous implementation
#copEV = convert(Lar.ChainOp, Lar.coboundary_0(EV))
#V = convert(Lar.Points, transpose(V))
#V, copEV, FE = Lar.Arrangement.planar_arrangement( V::Lar.Points, copEV::Lar.ChainOp )
	
	
	return V, copEV, FE
end


"""
	decomposition()::
	
Pairwise *intersection* in ``z = 0`` of *line segments* in ``σ ∪ I(σ)``, for each ``σ ∈ Sd−1``.

# Example 3D

```julia
V,FV,EV = model3d
model = model3d

```
"""
function decomposition(model::Lar.LAR)
	V,FV,EV = model
	dim = size(V,1)
	spatialindex = spaceindex(model)
	
	function submanifoldmap(vs)
		centroid = [sum(vs[k,:]) for k=1:size(vs,1)]/size(vs,2)
		# u1, u2 always independent
		u1 = normalize( centroid - vs[:,1] )
		u2 = normalize( vs[:,2] - vs[:,1] )
		u3 = normalize(cross(u1, u2))
		# u1, u2, u3 orthonormal
		u1 = cross(u2, u3)
		T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
		T[1:3,4] = - vs[:,1]
		R = Matrix{Float64}(LinearAlgebra.I, 4, 4)
		R[1:3, 1:3] = [u1 u2 u3]'
		return R*T  # roto-translation matrix
	end

	for Sigma in spatialindex
		sigma = Sigma[1]
		if dim == 3
			# transform Sigma s.t. Sigma[1], i.e. sigma, -> z=0
			vs = V[:, CV[sigma]]
			Q = submanifoldmap(vs)
			vq = Q * [vs; ones(1, size(vs,2))]
			v2d = vq[1:2,:]
	
		end
	end
	
end



"""
	Congruence()::
Graded bases of equivalence classes Ck (Uk ), with Uk = Xk /Rk for 0 ≤ k ≤ 2.

# Example

```julia
julia> 
```
"""
function Congruence 

end



"""
	Connection()::
Extraction of (X p , ∂p ), maximal connected components of Xd −1 (0 ≤ p ≤ h). d−1 d−1 +p

# Example

```julia
julia> 
```
"""
function Connection 

end



"""
	Bases()::
Computation of redundant cycle basis [∂d ] for each p-component, via TGW. 

# Example

```julia
julia> 
```
"""
function Bases 

end



"""
	Boundaries()::
Accumulation into H += [o]p (hole-set) of outer boundary cycle from each [∂d+]p . 

# Example

```julia
julia> 
```
"""
function Boundaries 

end



"""
	Containment()::
Computation of antisymmetric containment relation S between [o]p holes in H. 

# Example

```julia
julia> 
```
"""
function Containment 

end



"""
	Reduction()::
Transitive R reduction of S and generation of forest of flat trees ⟨[od ]p , [∂d ]p ⟩. 

# Example

```julia
julia> 
```
"""
function Reduction 

end



"""
	Adjoining()::
of roots [od ]r to (unique) outer cell, and non-roots [∂d+]q to container cells. 

# Example

```julia
julia> 
```
"""
function Adjoining 

end



"""
	Assembling()::
Quasi-block-diagonal assembly of matrices relatives to isolated components [∂d ]p . 

# Example

```julia
julia> 
```
"""
function Assembling 

end



"""
	Output()::
Global boundary map [∂d ] of A(Sd−1), and reconstruction of 0-chains of d-cells in Xd .

# Example

```julia
julia> 
```
"""
function Output 

end



