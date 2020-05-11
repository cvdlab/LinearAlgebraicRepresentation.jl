using Test
using SparseArrays

@testset "Model Generation" begin
    G = [0.0 0.0 1.0 1.0;  0.0 1.0 0.0 1.0]
    T1 = SparseArrays.sparse(Array{Int8,2}([1 1 0 0;  0 1 1 0;  0 0 1 1]))
    T2 = SparseArrays.sparse(Array{Int8,2}([1 1 1]))

    @test isnothing(CAGD.Model())

    m = CAGD.Model(G)
    @test size(m, 0, 1) == 2
    @test size(m, 0, 2) == 4
    @test size(m, 1, 1) == 0
    @test size(m, 1, 2) == 4
    @test size(m, 2, 1) == 0
    @test size(m, 2, 2) == 0
    @test_throws BoundsError size(m, 3, 1)

    m = CAGD.Model(G, [T1, T2])
    @test size(m, 0, 1) == 2
    @test size(m, 0, 2) == 4
    @test size(m, 1, 1) == 3
    @test size(m, 1, 2) == 4
    @test size(m, 2, 1) == 1
    @test size(m, 2, 2) == 3
    @test_throws BoundsError size(m, 3, 1)

    @test_throws ArgumentError CAGD.Model(G, [T1])
    @test_throws ArgumentError CAGD.Model(G, [T1, T1])
    @test_throws ArgumentError CAGD.Model(G, [T1, T1, T1])

end

