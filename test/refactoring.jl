using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation

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

data1 = (V,EV)
data2 = Lar.Struct([ Lar.t(2,2), Lar.r(pi/3), Lar.t(-1.5,-2.5), (W,EW) ])
data3 = Lar.Struct([ Lar.t(2,2), mycircle(2.5,16) ])
data4 = Lar.Struct([ Lar.t(3.5,3.5), mycircle(.25,16) ])
data5 = Lar.Struct([ Lar.t(5,3.5), mycircle(.5,16) ])
data6 = Lar.Struct([ Lar.t(5,3.5), mycircle(.25,16) ])

V,EV = input_collection( [ data1, data2, data3, data4, data5, data6 ] )

VV = [[k] for k in 1:size(V,2)];
Plasm.view( Plasm.numbering(.3)((V,[VV, EV])) )
```
Note that `V,EV` is not a cellular complex, since 1-cells intersect out of 0-cells.

# Example 3D

```julia
V,FV = Lar.sphere(radius=2)([3,4])
EV = Lar.simplexFacets(FV)
mysphere = V,FV,EV

data1 = mysphere
data2 = Lar.Struct([ Lar.t(0,1,0), mysphere ])
data3 = Lar.Struct([ Lar.t(0,0.5,0), Lar.s(0.4,0.4,0.4), mysphere ])
data4 = Lar.Struct([ Lar.t(4,0,0), Lar.s(0.8,0.8,0.8), mysphere ])
data5 = Lar.Struct([ Lar.t(4,0,0), Lar.s(0.4,0.4,0.4), mysphere ])
data = input_collection([ data1, data2, data3, data4, data5 ])
V,FV,EV = data

VV = [[k] for k in 1:size(data[1],2)];
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
	Indexing()::
Spatial index made by intersection of d interval-trees on bounding boxes of σ ∈ Sd −1 .

# Example

```julia
julia> 
```
"""
function Indexing 

end



"""
	Decomposition()::
Pairwise z = 0 intersection of line segments in σ ∪ I(σ), for each σ ∈ Sd−1.

# Example

```julia
julia> 
```
"""
function Decomposition 

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



