# [Lar.jl](@ref LinearAlgebraicRepresentation.Lar_header)

`Lar.jl` is a [Julia](http://julialang.org) library to perform geometrical computations on cellular complexes expressed through the [**Linear Algebraic Representation** (LAR)](https://www.sciencedirect.com/science/article/pii/S001044851300184X?via%3Dihub).
This library is developed and maintained by the [Computational Visual Design Laboratory (CVDLAB) of Università degli Studi di Roma Tre](https://github.com/cvdlab).

## Dependencies

`Lar.jl` has several Julia dependencies:

- [NearestNeighbors](https://github.com/KristofferC/NearestNeighbors.jl)
- DataStructures
- [IntervalTrees](https://github.com/BioJulia/IntervalTrees.jl)
- [Triangle](https://github.com/cvdlab/Triangle.jl)


## Docstrings conventions

**Bold** is used to point out theory concepts. For example, look at the 
"**1-skeletons**" word in the docstring of `Lar.skel_merge`:
```@docs
Lar.skel_merge(V1::Lar.Points, EV1::Lar.ChainOp, V2::Lar.Points, EV2::Lar.ChainOp)
```
`Monospace` is used for everything code related. Look this time at "`container`",
"`contained`" and "`Points`" in the docstring of `Lar.bbox_contains`:
```@docs
Lar.bbox_contains
```
!!! note
    In Julia REPL the `monospace` text is the one colored differently. In a terminal you will see something like:  
    ![Julia REPL monospace exaple](./images/monospace_juliarepl.png)
