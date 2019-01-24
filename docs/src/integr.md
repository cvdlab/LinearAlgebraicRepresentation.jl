# Finite integration of polynomials

A finite integration method, from "[Boundary integration over linear polyhedra](https://www.sciencedirect.com/science/article/pii/001044859090007Y)", is developed here, to compute various-order monomial integrals over polyhedral solids and surfaces in 3D space. The integration method can be used for the exact evaluation of domain integrals of trivariate polynomial forms.

## Integration Algorithms

Here we implement a finite solution
to both surface and volume integration of polynomials, by using a
triangulation of the domain boundary.   The evaluation of
surface and volume integrals is achieved by transformation into line
integrals over the boundary of every 2-simplex of a domain
triangulation.  A different approach to finite integration, using a
decomposition into volume elements induced by a boundary triangulation
is given in "[A symbolic method for calculating the integral properties of arbitrary nonconvex polyhedra](https://ieeexplore.ieee.org/document/6429334/)" where a closed formula for volume
integration over polyhedral volumes, by decomposing the solid into a set
of solid tetrahedra, but such a method cannot be used for surface
integrations.

Surface integrals are computed as a summation of integrals over a 
triangulation of the surface.  Any triangle is mapped into the unit 
triangle in the 2-space of parameters, where integrals 
of monomials become particularly 
simple.  Then formulae for integrals over polyhedral volumes are 
given.  They are easily derived by transforming volume integrals in 
surface integrals using the Divergence Theorem.  
It is possible to show that such integrals are 
computable in polynomial time, and that inertia moments are computable 
in $O(E)$ time, $E$ being the number of edges of the solid model of 
the integration domain.

An important feature of the integration formulae presented here is
that they can also be used with a partial model of a polyhedron,
consisting of the collection of its boundary's 2-cell cycles (loops).  
Loops are oriented
counter-clockwise if external, clockwise if internal to another loop. 
Such a  model, without explicit storage of face adjacencies, is
very frequently adopted in Computer Graphics.
In this case it is sufficient to consider any $(n+1)$-sided face (also
unconnected or multiply connected)  as topological sum of $n-1$
oriented triangles $t_i$, with vertices $\langle v_0, v_i,
v_{i+1}\rangle$, where $1\le i\le n-1$ .  In applying formulae
(\ref{19}) or (\ref{27}) to such a set of triangles, any edge that
does not belong to the original polygon will be computed twice, in the
two opposite directions.  These contributions to the whole integral
will mutually cancel each other, as they correspond to pairs of
line integrals evaluated along opposite paths

## Examples


### 2D integration

First some examples are given of the basic integration functions on 2D plane. The `II` integral, with parameters `alpha,beta,gamma = 0,0,0` returns the area of the domain `P`. It is important to notice that all domains (both 2D and 3D) are embedded in 3-space.

# Example  unit 3D triangle
```julia
julia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
3×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 0.0  0.0  0.0

julia> FV = [[1,2,3]]
1-element Array{Array{Int64,1},1}:
 [1, 2, 3]

julia> P = V,FV
([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])

julia> Lar.II(P, 0,0,0)
0.5
```
Then, more interesting examples are given. First, the surface of the affinely transformed (rotated and translated) unit square:

```julia
julia> V,FV = Lar.simplexGrid([1,1])
([0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0], Array{Int64,1}[[1, 2, 3], [2, 3, 4]])

julia> P = [V;[0 0 0 0]], FV
([0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0], 
Array{Int64,1}[[1, 2, 3], [2, 3, 4]])

julia> Lar.surface(P)
1.0

julia> p = Lar.Struct([Lar.t(0.5,0.5,0), Lar.r(0,0,pi/4), P]);

julia> q = Lar.struct2lar(p);

julia> Lar.surface(q)
1.0000000532124802
```
Then, the surface and the barycenter  of the simplicial grid `3 x 4` of unit 2-cells on the ``z=0`` plane is computed. Notice that the  `centroid` function cannot be used, since `P` is two-dimensional (flat) and embedded in 3-space. So, we use directly the centroid definition as first surface moments divided by area. Remember that the grid domain is ``3 x 4``.

```julia
julia> V,FV = Lar.simplexGrid([3,4]);

julia> P = [V; zeros(size(V,2))'], FV;

julia> Lar.surface(P)
12.0

julia> Lar.II(P,1,0,0)/Lar.II(P,0,0,0)
1.5

julia> Lar.II(P,0,1,0)/Lar.II(P,0,0,0)
2.0

julia> Lar.II(P,0,0,1)/Lar.II(P,0,0,0)
0.0
```


### 3D integration

The simplest example of volume integration is `volume` integral on the unit 3D tetrahedron `P`.

```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0

julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]
4-element Array{Array{Int64,1},1}:
 [1, 2, 4]
 [1, 3, 2]
 [4, 3, 1]
 [2, 3, 4]

julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], 
Array{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])

julia> Lar.volume(P)
0.16666666666666674
```

For more general polyhedrons, we must extract the boundary sourface, get a triangulation
and apply the 3D integration functions on the LAR of such models.


## Main Interface

```@docs
Lar.surface
```

```@docs
Lar.volume
```

```@docs
Lar.centroid
```
