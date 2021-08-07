# To be executed from CAGD.jl

todisplay = VERSION <= VersionNumber("1.2") ? true : false

using CAGD
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
if todisplay  using ViewerGL; GL = ViewerGL  end

include("./CAGD.jl/paper/examples/cellular_complex.jl");
include("./CAGD.jl/paper/examples/CSG_example.jl");
include("./CAGD.jl/paper/examples/cycles_boundaries.jl");
include("./CAGD.jl/paper/examples/planar_algebra_2D.jl");
include("./CAGD.jl/paper/examples/solid_algebra_3D.jl");
include("./CAGD.jl/paper/examples/space_arrangement_2d.jl");
include("./CAGD.jl/paper/examples/space_arrangement.jl");
include("./CAGD.jl/paper/examples/subdivision_of_2_cells.jl");
