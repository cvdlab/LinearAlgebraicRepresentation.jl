# LARLIB.jl

## Dependencies

`LARLIB.jl` has several Julia dependencies:

- [NearestNeighbors](https://github.com/KristofferC/NearestNeighbors.jl)
- DataStructures
- [IntervalTrees](https://github.com/BioJulia/IntervalTrees.jl)
- [TRIANGLE](https://github.com/furio/TRIANGLE.jl)


## Docstrings conventions

**Bold** is used to point out theory concepts. For example, look at the 
"**2-skeletons**" word in the docstring of `LARLIB.skel_merge`:
```@docs
LARLIB.skel_merge(V1::LARLIB.Points, EV1::LARLIB.Cells, V2::LARLIB.Points, EV2::LARLIB.Cells)
```
`Monospace` is used for everything code related. Look this time at "`container`",
"`contained`" and "`Points`" in the docstring of `LARLIB.bbox_contains`:
```@docs
LARLIB.bbox_contains
```
!!! note
    In Julia REPL the `monospace` text is the one colored differently. In a terminal you will see something like:  
    ![Julia REPL monospace exaple](./images/monospace_juliarepl.png)