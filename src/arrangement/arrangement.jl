module Arrangement
    using LARLIB
    using IntervalTrees
    using NearestNeighbors
    using TRIANGLE

    include("./minimal_cycles.jl")
    include("./dimension_travel.jl")
    include("./planar_arrangement.jl")
    include("./spatial_arrangement.jl")

    export planar_arrangement, spatial_arrangement

end