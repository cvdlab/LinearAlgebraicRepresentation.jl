# LARLIB.jl

`LARLIB.jl` is a [Julia](http://julialang.org) library to perform geometrical computations on cellular complexes expressed through the [**Linear Algebraic Representation** (LAR)](./lar.html).
This library is developed and maintained by the [Computational Visual Design Laboratory (CVDLAB) of Università degli Studi di Roma Tre](https://github.com/cvdlab).

## Dependencies

`LARLIB.jl` has several Julia dependencies:

- [NearestNeighbors](https://github.com/KristofferC/NearestNeighbors.jl)
- DataStructures
- [IntervalTrees](https://github.com/BioJulia/IntervalTrees.jl)
- [Triangle](https://github.com/cvdlab/Triangle.jl)


## Docstrings conventions

**Bold** is used to point out theory concepts. For example, look at the 
"**1-skeletons**" word in the docstring of `LARLIB.skel_merge`:
```@docs
LARLIB.skel_merge(V1::LARLIB.Points, EV1::LARLIB.ChainOp, V2::LARLIB.Points, EV2::LARLIB.ChainOp)
```
`Monospace` is used for everything code related. Look this time at "`container`",
"`contained`" and "`Points`" in the docstring of `LARLIB.bbox_contains`:
```@docs
LARLIB.bbox_contains
```
!!! note
    In Julia REPL the `monospace` text is the one colored differently. In a terminal you will see something like:  
    ![Julia REPL monospace exaple](./images/monospace_juliarepl.png)