# To be executed from CAGD.jl

todisplay = VERSION <= VersionNumber("1.2") ? true : false

using CAGD
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
if todisplay  using ViewerGL; GL = ViewerGL  end

include("cellular_complex.jl");
include("CSG_example.jl");
include("cycles_boundaries.jl");
include("planar_algebra_2D.jl");
include("solid_algebra_3D.jl");
include("space_arrangement_2d.jl");
include("space_arrangement.jl");
include("subdivision_of_2_cells.jl");
