using DataStructures


"""
	approxVal(PRECISION)(value)

Transform the float `value` to get a `PRECISION` number of significant digits.
"""
function approxVal(PRECISION)
    function approxVal0(value)
        out = round(value, digits=PRECISION)
        if out == -0.0
            out = 0.0
        end
        return out
    end
    return approxVal0
end





"""
	W,CW = simplifyCells(V,CV)

Find and remove the duplicated vertices and the incorrect cells.
Some vertices may appear two or more times, due to numerical errors
on mapped coordinates. Close vertices are identified, according to the
PRECISION number of significant digits.
"""
function simplifyCells(V,CV)
	PRECISION = 5
	vertDict = DefaultDict{Array{Float64,1}, Int64}(0)
	index = 0
	W = Array{Float64,1}[]
	FW = Array{Int64,1}[]

	for incell in CV
		outcell = Int64[]
		for v in incell
			vert = V[:,v]
			key = map(approxVal(PRECISION), vert)
			if vertDict[key]==0
				index += 1
				vertDict[key] = index
				push!(outcell, index)
				push!(W,key)
			else
				push!(outcell, vertDict[key])
			end
		end
		append!(FW, [[Set(outcell)...]])
	end
	return hcat(W...),FW
end



"""
	circle(radius=1.; angle=2*pi)(shape=36)

Compute an approximation of the circunference curve in 2D, centered on the origin.

With default values, i.e. `circle()()`, return the whole circonference of unit radius, approximated with a `shape=36` number of 1-cells.

# Example
```julia
julia> W,CW = Lar.circle()();

julia> GL.VIEW([
	GL.GLLines(W, CW, GL.COLORS[12]),
	GL.GLFrame
]);
```
"""
function circle(radius=1., angle=2*pi)
    function circle0(shape=[36])
        V, EV = cuboidGrid(shape)
        V = (angle/shape[1])*V
        V = hcat(map(u->[radius*cos(u); radius*sin(u)], V)...)
        W, EW = simplifyCells(V, EV)
        return W, EW
    end
    return circle0
end


"""
	helix(radius=1., pitch=1., nturns=2)(shape=36*nturns)

Compute the approximate elix curve in three-dimensional space, with basis on ``z=0`` plane and centered around the ``z`` axis. The `pitch` of a helix is the height of one complete helix `turn`, measured parallel to the axis of the helix.

# Example
```julia
julia> V, CV = Lar.helix(.1, .1, 10)()
([0.1 0.0984808 … 0.0984808 0.1; 0.0 0.0173648 … -0.0173648 0.0; 0.0 0.0027778 … 0.997222 1.0], Array{Int64,1}[[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11]  …  [351, 352], [352, 353], [353, 354], [354, 355], [355, 356], [356, 357], [357, 358], [358, 359], [359, 360], [360, 361]])

julia> GL.VIEW([
	GL.GLLines(V, CV, GL.COLORS[12]),
	GL.GLFrame
]);
```
"""
function helix(radius=1., pitch=1., nturns=2)
    function helix0(shape=36*nturns)
        angle = nturns*2*pi
        V, EV = cuboidGrid([shape])
        V = (angle/shape)*V
        V = hcat(map(u->[radius*cos(u);radius*sin(u);(pitch/(2*pi))*u], V)...)
        W, EW = simplifyCells(V, EV)
        return W, EW
    end
    return helix0
end


"""
	disk(radius=1., angle=2*pi)(shape=[36, 1])

Compute the cellular complex approximating a circular sector of 2D disk centered on the origin. In geometry, a disk is the region in a plane bounded by a circle. The `shape` array provides the number of approximating 2-cells.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.disk()()..., GL.COLORS[1],1 ),
	GL.GLFrame
]);
```
"""
function disk(radius=1., angle=2*pi)
    function disk0(shape=[36, 2])
        V, FV = simplexGrid(shape)
        V = [angle/shape[1] 0;0 radius/shape[2]]*V
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[v*cos(u);v*sin(u)] end, W)...)
        W, FW = simplifyCells(V, FV)
        FW = [cell for cell in FW if length(cell)==3]
		#EW = Lar.simplexFacets(FW)
        return W, FW #, EW
    end
    return disk0
end


"""
	helicoid(R=1., r=0.5, pitch=1., nturns=2)(shape=[36*nturns, 2])

Compute an approximation of the helicoid surface in 3D, with basis on ``z=0`` plane and centered around the ``z`` axis.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.helicoid()()..., GL.COLORS[1],1 ),
	GL.GLFrame
]);
```
"""
function helicoid(R=1., r=0.5, pitch=1., nturns=2)
    function helicoid0(shape=[36*nturns, 2])
        angle = nturns*2*pi
        V, CV = simplexGrid(shape)
        V = [angle/shape[1] 0;0 (R-r)/shape[2]]*V
        V = broadcast(+, V, [0, r])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[v*cos(u);v*sin(u);
        	(pitch/(2*pi))*u] end, W)...)
        return V, CV
    end
    return helicoid0
end


"""
	ring(r=1., R=2., angle=2*pi)(shape=[36, 1])

Compute the cellular 2-complex approximating a (possibly full) sector of a non-contractible disk. `R` and `r` are the external and the internal radiuses, respectively.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.ring()()..., GL.COLORS[1],1 ),
	GL.GLFrame
]);
```
"""
function ring(r=1., R=2., angle=2*pi)
    function ring0(shape=[36, 1])
		V, CV = cuboidGrid(shape)
		CV = [[[u,v,w],[w,v,t]] for (u,v,w,t) in CV]
		CV = reduce(append!,CV)
        V = [angle/shape[1] 0;0 (R-r)/shape[2]]*V
        V = broadcast(+, V, [0, r])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[v*cos(u);v*sin(u)] end, W)...)
        W, CW = simplifyCells(V, CV)
		CW = [cell for cell in CW if length(cell)==3]
		return W,CW
    end
    return ring0
end


"""
	cylinder(radius=.5, height=2., angle=2*pi)(shape=[36, 1])

Compute a cellular 2-complex, approximation of a right circular cylindrical surface in 3D. The open surface has basis on ``z=0`` plane and is centered around the ``z`` axis.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.cylinder()()..., GL.COLORS[1],1 ),
	GL.GLFrame
]);
```
"""
function cylinder(radius=.5, height=2., angle=2*pi)
    function cylinder0(shape=[36, 1])
        V, CV = Lar.cuboidGrid(shape)
		CV = [[[u,v,w],[w,v,t]] for (u,v,w,t) in CV]
		CV = reduce(append!,CV)
        V = [angle/shape[1] 0.0 ; 0.0 1.0/shape[2]]*V
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[radius*cos(u);radius*sin(u);
        	height*v] end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return cylinder0
end



"""
	sphere(radius=1., angle1=pi, angle2=2*pi)(shape=[18, 36])

Compute a cellular 2-complex, approximation of the two-dimensional closed surface, embedded in a three-dimensional Euclidean space. Geographical coordinates are user to compute the 0-cells of the complex.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.sphere()()..., GL.COLORS[1],0.75 ),
	GL.GLFrame
]);
```
"""
@enum surface triangled=1 single=2
function sphere(radius=1., angle1=pi, angle2=2*pi, surface=triangled)
    function sphere0(shape=[18, 36])
        V, CV = simplexGrid(shape)
        V = [angle1/shape[1] 0;0 angle2/shape[2]]*V
        V = broadcast(+, V, [-angle1/2, -angle2/2])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[radius*cos(u)*cos(v);
        	radius*cos(u)*sin(v);radius*sin(u)]end, W)...)
        W, CW = simplifyCells(V, CV)
        CW = [triangle for triangle in CW if length(triangle)==3]
        if Int(surface)==1
        	return W, CW
        elseif Int(surface)==2
        	return W,[collect(1:size(W, 2))]
        end
    end
    return sphere0
end



"""
	toroidal(r=1., R=2., angle1=2*pi, angle2=2*pi)(shape=[24, 36])

Compute a cellular 2-complex, approximation of the two-dimensional surface, embedded in a three-dimensional Euclidean space.
Toroidal is a closed surface having genus one, and therefore possessing a single "hole".
It can be constructed from a rectangle by gluing both pairs of opposite edges together with no twists.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.toroidal()()..., GL.COLORS[1],0.75 ),
	GL.GLFrame
]);
```
"""
function toroidal(r=1., R=2., angle1=2*pi, angle2=2*pi)
    function toroidal0(shape=[24, 36])
        V, CV = simplexGrid(shape)
        V = [angle1/(shape[1]) 0;0 angle2/(shape[2])]*V
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[(R+r*cos(u))*cos(v);
        	(R+r*cos(u))*sin(v);-r*sin(u)]end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return toroidal0
end



"""
	crown(r=1., R=2., angle=2*pi)(shape=[24, 36])

Compute a cellular 2-complex, approximation of a two-dimensional open surface, embedded in a three-dimensional Euclidean space.
This open surface is generated as an "half-torus", providing only the external shell.

# Example
```julia
julia> GL.VIEW([
	GL.GLGrid( Lar.crown()()..., GL.COLORS[1],0.75 ),
	GL.GLFrame
]);
```
"""
function crown(r=1., R=2., angle=2*pi)
    function crown0(shape=[12, 36])
        V, CV = simplexGrid(shape)
        V = [pi/shape[1] 0;0 angle/shape[2]]*V
        V = broadcast(+, V, [-pi/2, 0])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v)=p;[(R+r*cos(u))*cos(v);
       	(R+r*cos(u))*sin(v);-r*sin(u)]end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return crown0
end



"""
	cuboid(maxpoint::Array, full=false, minpoint::Array=zeros(length(maxpoint)))

Return a ``d``-dimensional cube, where ``d`` is the common length of arrays `minpoint` and
`maxpoint`.
If `flag=true` the cells of all dimensions (between 1 and ``d``) are generated.

```julia
julia> cuboid([-0.5, -0.5])
([0.0 0.0 -0.5 -0.5; 0.0 -0.5 0.0 -0.5], Array{Int64,1}[[1, 2, 3, 4]])

julia> cuboid([-0.5, -0.5, 0], true)
([0.0 0.0 … -0.5 -0.5; 0.0 0.0 … -0.5 -0.5; 0.0 0.0 … 0.0 0.0],
Array{Array{Int64,1},1}[Array{Int64,1}[[1], [2], [3], [4], [5], [6], [7], [8]],
Array{Int64,1}[[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], [6, 8], [1, 5], [2,
6], [3, 7], [4, 8]], Array{Int64,1}[[1, 2, 3, 4], [5, 6, 7, 8], [1, 2, 5, 6], [3, 4, 7,
8], [1, 3, 5, 7], [2, 4, 6, 8]], Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]]])

julia> V, (VV, EV, FV, CV) = Lar.cuboid([1,1,1], true);

julia> assembly = Lar.Struct([ (V, CV), Lar.t(1.5,0,0), (V, CV) ])

julia> GL.VIEW([
	GL.GLPol( Lar.struct2lar(assembly)..., GL.COLORS[1],0.75 ),
	GL.GLFrame ]);
```
"""
function cuboid(maxpoint::Array, full=false,
				minpoint::Array=zeros(length(maxpoint)))
	@assert( length(minpoint) == length(maxpoint) )
	dim = length(minpoint)
	shape = ones(Int, dim)
	cell = cuboidGrid(shape, full)
	size = maxpoint - minpoint
	out = apply(t(minpoint...) * s(size...), cell)
end



"""
	ball(radius=1, angle1=pi, angle2=2*pi)(shape=[18, 36,4])

Generate a cell decomposition of a *solid 3-sphere* in ``R^3``.
The variable `shape` provides the domain decomposition. Empty cells are removed after the *Cartesian -> Polar* coordinate mapping.
# Example
```julia
julia> GL.VIEW([
	GL.GLPol( Lar.ball()()..., GL.COLORS[1],0.5 ),
	GL.GLFrame ]);
```
"""
function ball(radius=1, angle1=pi, angle2=2*pi)
    function ball0(shape=[18, 36, 4])
        V, CV = cuboidGrid(shape)
        V = [angle1/shape[1] 0 0; 0 angle2/shape[2] 0; 0 0 radius/shape[3]]*V
        V = broadcast(+, V, [-(angle1)/2, -(angle2)/2, 0])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let
        (u, v, w)=p; [w*cos(u)*cos(v); w*cos(u)*sin(v); w*sin(u)] end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return ball0
end



"""
	rod(radius=1, height=3, angle=2*pi)(shape=[36, 1])

Compute a cellular 3-complex with a *single* 3-cell starting from a cyclindrical surface generated with the same parameters.

# Example
```julia
julia> rod()()[1]
# output
3×74 Array{Float64, 2}:
 -0.34202   0.984808   1.0          1.0  …   0.984808  -0.866025  -1.0           0.766044
  0.939693  0.173648  -2.44929e-16  0.0     -0.173648  -0.5        1.22465e-16  -0.642788
  3.0       3.0        0.0          0.0      0.0        3.0        3.0           0.0

julia> rod()()[2]
# output
1-element Array{Array{Int64, 1}, 1}:
 [0, 1, 2, 3, 4, 5, 6, 7, 8, 9  …  64, 65, 66, 67, 68, 69, 70, 71, 72, 73]

julia> GL.VIEW([
 	GL.GLPol( Lar.rod()()..., GL.COLORS[1],0.5 ),
 	GL.GLFrame ]);
```
"""
function rod(radius=1., height=3., angle=2*pi)
    function rod0(shape=[36, 1])
        V, CV = cylinder(radius, height, angle)(shape)
        return V, [collect(1:size(V, 2))]
    end
    return rod0
end



"""
	hollowCyl(r=1., R=2., height=6., angle=2*pi)(shape=[36, 1, 1])

Compute the cellular 3-complex approximating a solid cylinder with a
internal axial hole. The model is meshed with cubical 3-cells.

# Example
```julia
julia> GL.VIEW([
 	GL.GLPol( Lar.hollowCyl()()..., GL.COLORS[1],0.5 ),
 	GL.GLFrame ]);
```
"""
function hollowCyl(r=1., R=2., height=6., angle=2*pi)
    function hollowCyl0(shape=[36, 1, 1])
        V, CV = cuboidGrid(shape)
        V = [angle/shape[1] 0 0;0 (R-r)/shape[2] 0;0 0 height/shape[3]]*V
        V = broadcast(+, V, [0, r, 0])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v, z)=p;[v*cos(u);v*sin(u);z] end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return hollowCyl0
end



"""
	hollowBall(r=1., R=2., angle1=pi, angle2=2*pi)(shape=[36, 1, 1])

Compute the cellular 3-complex approximating a 3-sphere. The model is meshed with cubical 3-cells, where the mesh has default decomposition size `[24, 36, 8]`.

# Example
```julia
julia> V, CV = Lar.hollowBall(1, 2, pi/2, pi/2)([6, 12, 4]);

julia> GL.VIEW([
 	GL.GLPol( V, CV, GL.COLORS[1],0.5 ),
 	GL.GLFrame ]);
...
```

"""
function hollowBall(r=1., R=1., angle1=pi, angle2=2*pi)
    function hollowBall0(shape=[24, 36, 3])
        V, CV = cuboidGrid(shape)
        V = [angle1/shape[1] 0 0; 0 angle2/shape[2] 0; 0 0 (R-r)/shape[3]]*V
        V = broadcast(+, V, [-(angle1)/2, -(angle2)/2, r])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let
        (u, v, w)=p; [w*cos(u)*cos(v); w*cos(u)*sin(v); w*sin(u)] end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return hollowBall0
end



"""
	torus(r=1., R=2., h=.5, angle1=2*pi, angle2=2*pi)(shape=[24, 36, 4])

Compute the cellular 3-complex approximating the solid torus in 3D. The model is meshed with cubical 3-cells, where the mesh has default decomposition size `[24, 36, 4]`. See also: [`toroidal`](@toroidal). `h` is radius of the circular hole inside the solid.

# Example
```julia
julia> GL.VIEW([
 	GL.GLPol( Lar.torus(1., 2., .5, pi, pi)()..., GL.COLORS[1],0.5 ),
 	GL.GLFrame ]);
```
"""
function torus(r=1., R=2., h=.5, angle1=2*pi, angle2=2*pi)
    function torus0(shape=[24, 36, 4])
        V, CV = cuboidGrid(shape)
        V = [angle1/shape[1] 0 0;0 angle2/shape[2] 0;0 0 r/shape[3]]*V
        V = broadcast(+, V, [0, 0, h])
        W = [V[:, k] for k=1:size(V, 2)]
        V = hcat(map(p->let(u, v, z)=p;[(R+z*cos(u))*cos(v);(R+z*cos(u))*sin(v);
        	-z*sin(u)] end, W)...)
        W, CW = simplifyCells(V, CV)
        return W, CW
    end
    return torus0
end



#"""
#	pizza(r=.1, R=1., angle=pi)(shape=[24, 36])
#
#Compute a cellular 3-complex with a single convex 3-cell.
#
## Example
#```julia
#julia> model = pizza(r=.1, R=1., angle=pi)([12,18])
#
#julia> model[1]
#3×249 Array{Float64,2}:
# 0.799215  0.997552   0.871665   0.573303  …   0.8463      0.964127  0.985637    0.0   0.0
# 0.670621  0.175895   0.573303   0.871665      0.55662     0.415884  0.2336      0.0   0.0
# 0.025     0.0482963  0.025     -0.025        -0.0482963  -0.0       0.0482963  -0.05  0.05
#
#julia> model[2]
#1-element Array{Array{Int64,1},1}:
# [1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  240, 241, 242, 243, 244, 245, 246, 247, 248, 249]
#
#julia> using Plasm
#
#julia> Plasm.view(model)
#```
#"""
#function pizza(r=.1, R=1., angle=pi)
#    function pizza0(shape=[24, 36])
#        V, CV = crown(r, R, angle)(shape)
#        W = [Any[V[h, k] for h=1:size(V, 1)] for k=1:size(V, 2)]
#        X = hcat(collect(Set(W))...)
#        V = hcat(X, [0 0;0 0;-r r])
#        return V, [collect(1:size(V, 2))]
#    end
#    return pizza0
#end
