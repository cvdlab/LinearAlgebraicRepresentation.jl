using Base.Test
using LinearAlgebraicRepresentation


@testset "submanifold_mapping test" begin
    V = rand(3, 3)
    m = LinearAlgebraicRepresentation.Arrangement.submanifold_mapping(V)
    err = LinearAlgebraicRepresentation.ERR 
    
    @test any(map((x,y)->-err<x-y<err, m*inv(m), eye(4)))
    @test any(x->-err<x<err, ([V [1; 1; 1]]*m)[:, 3])
end

