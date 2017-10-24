module LARLIB

    try
        Pkg.installed("TRIANGLE")
    catch
        Pkg.clone("https://github.com/furio/TRIANGLE.jl.git")
        Pkg.build("TRIANGLE")
    end

    module Arrangement
        const Verts = Array{Float64, 2}
        const Cells = SparseMatrixCSC{Int8, Int}
        const Cell = SparseVector{Int8, Int}
        
        include("./planar_arrangement.jl")
        include("./spatial_arrangement.jl")
    end

    function skel_merge(V1, EV1, V2, EV2)
        Arrangement.skel_merge(V1, EV1, V2, EV2)
    end

    function skel_merge(V1, EV1, FE1, V2, EV2, FE2)
        Arrangement.skel_merge(V1, EV1, FE1, V2, EV2, FE2)
    end

    function spatial_arrangement(V, EV, FE)
        Arrangement.spatial_arrangement(V, EV, FE)
    end

    function planar_arrangement(V, EV)
        Arrangement.planar_arrangement(V, EV)
    end
end
