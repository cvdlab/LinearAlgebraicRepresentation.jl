const Verts = Array{Float64, 2}
const Cells = SparseMatrixCSC{Int8, Int}
const Cell = SparseVector{Int8, Int}

include("./planar_arrangement.jl")

include("./dimension_travel.jl")
include("./utilities.jl")
