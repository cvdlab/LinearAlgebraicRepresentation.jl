module LARLIB
    module PlanarArrangement
        include("./planar_arrangement.jl")
    end

    function planar_arrangement(V, EV)
        PlanarArrangement.planar_arrangement(V, EV)
    end
end
