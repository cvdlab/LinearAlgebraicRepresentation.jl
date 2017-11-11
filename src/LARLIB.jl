module LARLIB

    using NearestNeighbors
    using NearestNeighbors
    using IntervalTrees
    using TRIANGLE
    
    const Verts = Array{Float64, 2}
    const Cells = SparseMatrixCSC{Int8, Int}
    const Cell = SparseVector{Int8, Int}
    

    include("./utilities.jl")
    include("./minimal_cycles.jl")
    include("./dimension_travel.jl")
    include("./planar_arrangement.jl")
    include("./spatial_arrangement.jl")
    
end
