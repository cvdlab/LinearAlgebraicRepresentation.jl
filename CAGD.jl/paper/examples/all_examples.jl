todisplay = true

using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
using SparseArrays
if todisplay  using ViewerGL; GL = ViewerGL  end

include("./updatedLarFunctions.jl");


include("./cellular_complex.jl");
include("./CSG_example.jl");
include("./cycles_boundaries.jl");
include("./planar_algebra_2D.jl");
include("./solid_algebra_3D.jl");
include("./space_arrangement_2d.jl");
include("./space_arrangement.jl");
include("./subdivision_of_2_cells.jl");
