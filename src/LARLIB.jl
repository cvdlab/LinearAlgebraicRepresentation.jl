module LARLIB

    using NearestNeighbors
    using DataStructures
    using NearestNeighbors
    using IntervalTrees
    using TRIANGLE
    
    const Verts = Array{Float64, 2}
    const Cells = SparseMatrixCSC{Int8, Int}
    const Cell = SparseVector{Int8, Int}
    const LarCells = Array{Array{Int, 1}, 1}
    

    include("./utilities.jl")
    include("./minimal_cycles.jl")
    include("./dimension_travel.jl")
    include("./planar_arrangement.jl")
    include("./spatial_arrangement.jl")
    include("./largrid.jl")
    
end
