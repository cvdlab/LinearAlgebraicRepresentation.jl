const Verts = Array{Float64, 2}
const Cells = SparseMatrixCSC{Int8, Int}
const Cell = SparseVector{Int8, Int}

try
   Pkg.installed("TRIANGLE")
catch
   Pkg.clone("https://github.com/furio/TRIANGLE.jl.git")
   Pkg.build("TRIANGLE")
end

include("./planar_arrangement.jl")
include("./dimension_travel.jl")
include("./utilities.jl")
