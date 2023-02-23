using SparseArrays
using ViewerGL
GL = ViewerGL
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

function allexamples2d()
   include("twosquares.jl")
   include("arrangement2d.jl")
   include("biconnectedgraphs.jl")
   include("graphalgorithms.jl")
   include("orient2d.jl")
   include("point_in_polygon.jl")
   include("randomarrangement2d.jl")
   include("randombubbles.jl")
   include("randomlines.jl")
   include("randomshapes.jl")
   include("stresstest2d.jl")
   include("svg2lar.jl")
end

allexamples2d()
