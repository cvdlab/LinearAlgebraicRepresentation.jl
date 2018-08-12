using Test
using LARLIB

@testset "submanifold_mapping test" begin
    V = rand(3, 3)
    m = LARLIB.Arrangement.submanifold_mapping(V)
    err = 1e-10

    @test any(map((x,y)->-err<x-y<err, m*inv(m), Matrix{Float64}(I, 4, 4)))
    @test any(x->-err<x<err, ([V [1; 1; 1]]*m)[:, 3])
end

