# Interface to cell and chain complexes

Most part of text in this page is derived from []() and from []().

## Introduction

With increased complexity of geometric data, **topological models** play an increasingly important role beyond *boundary representations*, *assemblies*, *finite elements*, *image processing*, and other traditional modeling applications. While many graph- and index- based data structures have been proposed, no standard representation has emerged as of now. Furthermore, such representations typically do not deal with representations of mappings and functions and do not scale to support parallel processing, open source, and client-based architectures. 

A proper mathematical model for all topological structures is a **(co)chain complex**: a sequence of linear **(co)chain spaces** and linear **(co)boundary mappings**. This in turn implies all topological structures may be represented by **a collection of sparse matrices**. We propose a **Linear Algebraic Representation (LAR)** scheme for ``mod 2`` (co)chain complexes using CSC sparse matrices and show that it supports variety of topological computations using standard matrix algebra, without any overhead in space or running time.
With the LAR scheme, only the *characteristic functions* (see below) of ``d``-cells as vertex subsets are necessary for representing *polytopal complexes*. Examples include simplicial, cubical, and Voronoi complexes.

## Data structures

All types of cell complexes and functions over cell complexes are properly represented by a (co)chain complex, that captures all combinatorial relationships of interest in solid and physical modeling formally and unambiguously. According to classical results from algebraic topology techniques, a (co)chain complex and all associated combinatorial operations are readily represented using *standard techniques* from linear algebra, giving rise to a *Linear Algebraic Represention* (LAR) scheme.

In this package, we provide LAR data structures and algorithms using *compressed sparse column (CSC) matrices*, that introduce no computational overhead and are asymptotically as efficient as (and usually better than) many other popular topological data structures. Our aim is to provide a representation that supports all topological constructions and queries that arise in typical *cellular decomposition* of space (mesh, image, boundary, etc).

An **arrangement** is the *decomposition of d-dimensional space* into connected and *relatively open cells* of lower dimensions, induced by an intersection of a finite collection of geometric objects. A *planar collection* S may include line segments, open or closed polygonal lines, polygons, two-dimensional meshes, and discrete images in 2D. A *space collection* may include 3D polygons, polygonal meshes, B-reps of solid models—either manifold or non-manifold, three-dimensional CAE meshes, and volumetric images in 3D.

In this package, we have implemented the *computation of the arrangement* produced by a  *set of cellular complexes* in either 2D or 3D. Our goal is to provide a complete description of the plane or space decomposition induced by the input, into cells of dimensions 0, 1, 2 or 3.

### Characteristic matrices

A precise mathematical definition of a cellular complex is not trivial; we may rely on the intuitive idea of constructing a space *by gluing together* a number of *building blocks* of different dimensions, called **cells**.

The **characteristic function** ``\chi _A : S \to \{ 0, 1 \}`` is a function defined on a set ``S = \{s_j\}``, that indicates membership of an element ``s_j`` in a subset ``A \subseteq S``, having the value 1 for all elements of ``A`` and the value 0 for all elements of ``S`` not in ``A``.
We call **characteristic matrix** ``M`` of a collection of subsets ``A_i \subseteq S``  ``(i=1,...,n)`` the binary matrix ``M=(m_{ij})``, with ``m_{ij} = \chi_{A_i}(s_j)``.

#### Examples

Binary matrix representing by rows the `p`-cells of a cellular complex.
The input parameter must be of `Cells` type. Return a sparse binary matrix, 
providing the basis of a ``Chain`` space of given dimension. Notice that the 
number of columns is equal to the number of vertices (0-cells). 

First the cellular complex describing the 0-, 1-, 2-, anf 3-faces of a single
unit cube is generated. 

```julia
V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true); 
```

Then, we may see the characteristic matrix of 1-cells (edges), with two ones per row:

```julia
julia> Lar = LinearAlgebraicRepresentation

julia> Matrix(Lar.characteristicMatrix(EV))
12×8 Array{Int8,2}:
 1  1  0  0  0  0  0  0
 0  0  1  1  0  0  0  0
 0  0  0  0  1  1  0  0
 0  0  0  0  0  0  1  1
 1  0  1  0  0  0  0  0
 0  1  0  1  0  0  0  0
 0  0  0  0  1  0  1  0
 0  0  0  0  0  1  0  1
 1  0  0  0  1  0  0  0
 0  1  0  0  0  1  0  0
 0  0  1  0  0  0  1  0
 0  0  0  1  0  0  0  1
```
and of 2-cells (faces):

```julia
julia> Matrix(Lar.characteristicMatrix(FV))
6×8 Array{Int8,2}:
 1  1  1  1  0  0  0  0
 0  0  0  0  1  1  1  1
 1  1  0  0  1  1  0  0
 0  0  1  1  0  0  1  1
 1  0  1  0  1  0  1  0
 0  1  0  1  0  1  0  1
```

Finally, the boundary of the single 3-cell contains all the 0-cells (vertices). Of course, the 3D cube has 12 edges in `EV`, 6 faces in `FV`, and one 3-cell in `CV`:

```julia
julia> Matrix(Lar.characteristicMatrix(CV))
1×8 Array{Int8,2}:
 1  1  1  1  1  1  1  1
```



### Chain bases

In algebraic topology, a ``k``-chain is a *formal linear combination* of the ``k``-cells in a cell complex. In *simplicial* complexes (respectively, *cubical* complexes), ``k``-chains are combinations of ``k``-simplices (respectively, ``k``-cubes).

Let ``\sigma`` be an oriented cell in ``X`` and ``\mu \in G``. The elementary chain whose value is ``\mu`` on ``\sigma``, ``-\mu`` on ``-\sigma`` and ``0`` on any other cell in ``X`` is denoted ``\mu\sigma`` . Each chain can then be written in a *unique way* as a *sum of elementary chains*. With abuse of notation, we do NOT distinguish between *cells* and *singleton chains* (i.e., the elementary chains whose value is ``1\sigma`` for some cell ``\sigma``), used as elements of the **standard bases** of chain groups.

Chains are often thought of as *attaching orientation* and *multiplicity* to cells: if coefficients are extracted from the group ``G = (\{ −1, 0, 1 \}, +) ≃ (\mathbf{Z}_3, +)``, then cells can only be discarded or selected, possibly inverting their orientation. 
A *``p``-cycle* is a *closed ``p``-chain*, i.e. a ``p``-chain *without boundary*. 
It is useful to select a conventional choice to orient the singleton chains (*single cells*) automatically. 0-cells are considered all positive. The ``p``-cells, for ``1 ≤ p ≤ d-1``, can be given an *internal orientation* according to the orientation of the first (``p − 1``)-cell in their *canonical representation*, i.e. sorted on indices of their (``p − 1``)-cycle. Finally, a ``d``-cell may be oriented as the sign of its ``oriented volume``.


#### Examples

A compact representation of bases of ``p``-cells is provided by `Cells` type, defined as `Array{Array{Int,1}}`, where each element codifies a cell as the array of indices to vertices on the boundary of the cell:

```julia
julia> V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true);

julia> FV
6-element Array{Array{Int64,1},1}:
 [1, 2, 3, 4]
 [5, 6, 7, 8]
 [1, 2, 5, 6]
 [3, 4, 7, 8]
 [1, 3, 5, 7]
 [2, 4, 6, 8]
```

A simplicial decomposition of the unit cube with six 3-cells (tetrahedra), and a simplicial decomposition of the domain ``[0,5] \times [0,1]`` with ten 2-cells (triangles) follows:

```julia
julia> V,CV = Lar.simplexGrid([1,1,1]);

julia> V
3×8 Array{Float64,2}:
 0.0  1.0  0.0  1.0  0.0  1.0  0.0  1.0
 0.0  0.0  1.0  1.0  0.0  0.0  1.0  1.0
 0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0

julia> CV   # bases of tetrahedra
6-element Array{Array{Int64,1},1}:
 [1, 2, 3, 5]
 [2, 3, 5, 6]
 [3, 5, 6, 7]
 [2, 3, 4, 6]
 [3, 4, 6, 7]
 [4, 6, 7, 8]

julia> W,FW = Lar.simplexGrid([5,1]);

julia> W
2×12 Array{Float64,2}:
 0.0  1.0  2.0  3.0  4.0  5.0  0.0  1.0  2.0  3.0  4.0  5.0
 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> FW   # bases of triangles
10-element Array{Array{Int64,1},1}:
 [1, 2, 7]  
 [2, 7, 8]  
 [2, 3, 8]  
 [3, 8, 9]  
 [3, 4, 9]  
 [4, 9, 10] 
 [4, 5, 10] 
 [5, 10, 11]
 [5, 6, 11] 
 [6, 11, 12]
```
It is worthwhile to note that the above examples provide ``p``-bases of suitable dimensions, equal to the number of ``p``-cells,  for the corresponding chain complexes.   


### (Co)boundary operators

**Boundary operators** are maps ``\partial_p : C_p \to C_{p−1}`` between chain spaces, i.e. between spaces of subsets of cells with different dimension, with ``1 \leq p \leq d``, hence for a *cellular 2-complex* we have two operators, denoted as ``\partial_2 : C_2 \to C_1`` and ``\partial_1 : C_1 \to C_0``, respectively. Since they are linear maps between linear spaces, may be represented by matrices of coefficients ``[\partial_2]`` and ``[\partial_1]`` from the corresponding groups. We use the groups ``\{0, 1\}`` and ``\{-1, 0, 1\}`` for *unsigned* and *signed* coefficients, respectively.

The concept of *cochain*  in a group ``C^p`` of linear maps from chains ``C_p \to \Re`` allows for the *association of numbers* not only to *single cells*, as done by chains, but also to *assemblies of cells*. A cochain is hence the association of every discretized subdomain (chain) of a cell complex with a numeric quantity, usually resulting from a *discrete integration* over a chain.

**Coboundary operators** are maps ``\delta^p : C^{p} \to C^{p+1}``, with each linear space  ``C^p`` of ``p``-cochains isomorphic to the space ``C_p`` of ``p``-chain. Therefore in this package only Chain spaces are used. Notice that ``[\partial_p] = [\delta^p]^t``.
This property is often utilized in our algorithms.


#### Examples

```julia
julia> V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true);

julia> EV
12-element Array{Array{Int64,1},1}:
 [1, 2]
 [3, 4]
   ...
 [3, 7]
 [4, 8]

julia> Lar.boundary_1( EV::Lar.Cells )
8×12 SparseMatrixCSC{Int8,Int64} with 24 stored entries:
  [1 ,  1]  =  -1
  [2 ,  1]  =  1
  [3 ,  2]  =  -1
	...       ...
  [7 , 11]  =  1
  [4 , 12]  =  -1
  [8 , 12]  =  1

julia> Matrix(Lar.boundary_1(EV::Cells))
8×12 Array{Int8,2}:
 -1   0   0   0  -1   0   0   0  -1   0   0   0
  1   0   0   0   0  -1   0   0   0  -1   0   0
  0  -1   0   0   1   0   0   0   0   0  -1   0
  0   1   0   0   0   1   0   0   0   0   0  -1
  0   0  -1   0   0   0  -1   0   1   0   0   0
  0   0   1   0   0   0   0  -1   0   1   0   0
  0   0   0  -1   0   0   1   0   0   0   1   0
  0   0   0   1   0   0   0   1   0   0   0   1
```
Notice that the matrix ``[\partial_1]``,  generated by the function `boundary_1` applied to the 1-cell basis `EV`, contains two non-zero elements per column, where the associated edge (1-cell) is oriented from the vertex (row) of lesser index, towards the vertex of greater index, according to our numbering convention.


###  Chain complexes

A **chain complex**, for our purposes, is an algebraic structure that consists of *a sequence of linear spaces*  and a *sequence of linear maps* between consecutive linear spaces,  such that the **image** of each map (subspace of *boundaries* of ``p``-chains) is included in the **kernel** (subspace of *cycles* of ``(p-1)``-chains) of the next. 

The set of all ``k``-chains forms a group and the *sequence of these groups* is called a **chain complex**. In computing the arrangement ``A(S)`` induced by ``S``, we actually compute the whole chain complex ``C\bullet`` generated by the cell complex ``X := A(S)``. For example, in 3D we compute all objects and arrows (morphisms) in the diagram below, and hence we obtain a computational knowledge of space subdivision homology, including the Euler number. 
 

#### Examples

From the minimal possible input, construct the whole
two-dimensional chain complex, i.e. the bases for linear spaces C_1 and 
C_2, of 1-chains and  2-chains, and the signed coboundary operators from 
C_0 to C_1 and from C_1 to C_2.

##### 2D Chain complex

Start with the 1-skeleton (set of 1-cells) of a 2D small cuboidal grid (made of squares);
in other words, suppose we only know the edges of the grid:

```julia
julia> W = 
 [0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0  3.0  3.0  3.0  3.0
  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0]
# output  
 2×16 Array{Float64,2}: ...

julia> EW = 
[[1, 2],[2, 3],[3, 4],[5, 6],[6, 7],[7, 8],[9, 10],[10, 11],[11, 12],[13, 14],
 [14, 15],[15, 16],[1, 5],[2, 6],[3, 7],[4, 8],[5, 9],[6, 10],[7, 11],[8, 12],
 [9, 13],[10, 14],[11, 15],[12, 16]]
# output  
24-element Array{Array{Int64,1},1}: ...
```
We go to compute the arrangement of the 2D space induced by the above, i.e. the full **chain complex** generated by `(W,EW)`. This one is returned by the evaluation of the expression `chaincomplex(W,EW)`. The output variable `bases`  will contain the meaningful cell bases, i.e. those of dimension 1 and 2, since dimension 0 -- isolated 0-cells -- is not so.

```julia
julia> V,bases,coboundaries = Lar.chaincomplex(W,EW)

julia> bases[1]	# edges
24-element Array{Array{Int64,1},1}: ...

julia> bases[2] # faces -- previously unknown !!
9-element Array{Array{Int64,1},1}: ...
```
Analogously, the `coboundaries` variable will contain the ``[\delta_1]`` and ``[\delta_2]`` matrices, of type `SparseMatrixCSC{Int8,Int64}`

```julia
julia> coboundaries[1] # coboundary_1 
24×16 SparseMatrixCSC{Int8,Int64} with 48 stored entries: ...

julia> Matrix(coboundaries[2]) # coboundary_1: faces as oriented 1-cycles of edges
9×24 Array{Int8,2}:
 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0
  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0
  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0
  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0
  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0
  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0
  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0
  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0  0
  0  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1
```
Note that the last matrix contains by rows the 2-cycles corresponding to the (previously) unknown 2-basis `FV` that can now be easily computed.
Notice also that columns corresponding to *interior edges* (1-cells) contain two non-zeros *of opposite sign*. Hence the computed 2-complex is **coherently oriented** by the matrix rows, actually by construction. 

The 2-boundary operator matrix, i.e. `transpose(coboundaries[2])`, can be used to compute the boundary of every possible 2-chain, by matrix multiplication times the coordinate (binary) representation of the 2-chain, implemented by the type `Chain`, defined as `SparseVector{Int8, Int}`.



##### 3D Chain complex

The example discussed here concerns two unit cubes in 3D, where the second is rotated and translated, up to intersect only partially the firat cube. First we prepare our data, using a very simple hierarchical aggregation via a `Struct` object, to get a representation of faces and edges of both cubes in two `Cells` arrays `FV` and `EV`.

```julia
julia> cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )

julia> cube_2 = Lar.Struct([Lar.t(0,0,0.5), Lar.r(0,0,pi/3), cube_1])

julia> V,FV,EV = Lar.struct2lar(Lar.Struct([ cube_1, cube_2 ]))
```
Then we compute the 3D space arrangement induced by `FV`, providing the auxiliary information in `EV`, and getting back `V,`bases`,coboundaries`.
Both `bases` and `coboundaries` are then disassembled into their component data structures. The *actual discoveries* computed by the `arrangement` algorithms, called by the `chaincomplex` function, are the NEW `EV`, `FV`, `CV` basis and the `cscEV, cscFE, cscCF` operators, which stands  for ``[\delta_0]``,``[\delta_1]``, and``[\delta_2]``, i.e. the **computations of the solid 3D cells** generated by the arrangement of space, including their **full topology**.


```julia
julia> V,bases,coboundaries = Lar.chaincomplex(V,FV,EV)

julia> (EV, FV, CV), (cscEV, cscFE, cscCF) = bases,coboundaries

julia> FV # bases[2]
18-element Array{Array{Int64,1},1}:
 [1, 3, 4, 6]            
 [2, 3, 5, 6]            
 [7, 8, 9, 10]           
 [1, 2, 3, 7, 8]         
 [4, 6, 9, 10, 11, 12]   
 [5, 6, 11, 12]          
 [1, 4, 7, 9]            
 [2, 5, 11, 13]          
 [2, 8, 10, 11, 13]      
 [2, 3, 14, 15, 16]      
 [11, 12, 13, 17]        
 [11, 12, 13, 18, 19, 20]
 [2, 3, 13, 17]          
 [2, 13, 14, 18]         
 [15, 16, 19, 20]        
 [3, 6, 12, 15, 19]      
 [3, 6, 12, 17]          
 [14, 16, 18, 20]        

julia> CV # bases[3]
3-element Array{Array{Int64,1},1}:
 [2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20]
 [2, 3, 5, 6, 11, 12, 13, 17]                    
 [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]    
```	
Let note that the LAR of `FV` includes faces with 5 and 6 vertices, even non convex.  Also, the `CV` variable contains the LAR of the three solid parts the two cubes are split into. 

```julia
julia> cscEV # coboundaries[1]
34×20 SparseMatrixCSC{Int8,Int64} with 68 stored entries: ...

julia> cscFE # coboundaries[2]
18×34 SparseMatrixCSC{Int8,Int64} with 80 stored entries: ...

julia> cscCF # coboundaries[3]
4×18 SparseMatrixCSC{Int8,Int64} with 36 stored entries: ...
```	


## Main Interface

```@docs
Lar.characteristicMatrix
```

```@docs
Lar.boundary_1
```

```@docs
Lar.coboundary_0
```

```@docs
Lar.u_coboundary_1
```

```@docs
Lar.u_boundary_2
```

```@docs
Lar.coboundary_1
```

```@docs
Lar.chaincomplex
```
