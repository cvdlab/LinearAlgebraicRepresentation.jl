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

The **characteristic function** ``\chi_A : S → \{0, \1}`` is a function defined on a set ``S = \{s_j\}``, that indicates membership of an element ``s_j`` in a subset ``A ⊆ S``, having the value 1 for all elements of ``A`` and the value 0 for all elements of ``S`` not in ``A``.
We call **characteristic matrix** ``M`` of a collection of subsets ``Ai ⊆ S (i=1,...,n)`` the binary matrix ``M=(m_{ij})``, with ``m_{ij} = \chi_{A_i}(s_j)``.

#### Examples




### Chain bases

Let σ be an oriented cell in X and д ∈ G. The elementary chain whose value is д on σ, −д on −σ and 0 on any other cell in X is denoted дσ . Each chain can then be written in a unique way as a sum of elementary chains. With abuse of notation, we do not distinguish between cells and singleton chains (i.e., the elementary chains whose value is 1 σ for some cell σ ), used as elements
of the standard bases7 of chain groups.

Chains are often thought of as attaching orientation and multiplicity to cells: if coefficients are
taken from the group G = ({−1, 0, 1}, +) ≃ (Z3, +), then cells can only be discarded or selected, possibly inverting their orientation (see [20]). A p-cycle is a closed p-chain, i.e. a p-chain without boundary. It is useful to select a conventional choice to orient the singleton chains (single cells) automatically. 0-cells are considered all positive. The p-cells, for 1 ≤ p ≤ d − 1, can be given a coherent (internal) orientation according to the orientation of the first (p − 1)-cell in their canonical representation sorted on indices of their (p − 1)-cycle. Finally, a d-cell may be oriented as the sign of its oriented volume.


#### Examples




### (Co)boundary operators

The boundary operators are maps Cp → Cp−1, with 1 ≤ p ≤ d, hence for a 2-complex we have two operators, denoted as ∂2 : C2 → C1 and ∂1 : C1 → C0, respectively. Since they are linear maps between linear spaces, may be represented by matrices of coefficients [∂2] and [∂1] from the corresponding groups. 

The concept of cochain φp in a group Cp of linear maps from chains Cp to R allows for the association of numbers not only to single cells, as done by chains, but also to assemblies of cells. A cochain is hence the association of discretized subdomains of a cell complex with a global numeric quantity, usually resulting from a discrete integration over a chain.


#### Examples




###  Chain complexes

In mathematics, a chain complex is an algebraic structure that consists of a sequence of abelian groups (or modules) and a sequence of homomorphisms between consecutive groups such that the image of each homomorphism is included in the kernel of the next. 

 in computing the arrangement A(S) induced by S, we actually compute the whole chain complex C• generated by the cell complex X := A(S). For example, in 3D we compute all objects and arrows (morphisms) in the diagram below, and hence we obtain a computational knowledge of space subdivision homology [36, 37], including the Euler number. 
 

#### Examples





