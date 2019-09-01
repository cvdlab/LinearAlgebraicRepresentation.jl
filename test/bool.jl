using Test
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

@testset "bool.jl file Tests" begin

    n,m,p = 1,1,1
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
    cube = V,FV,EV

    threecubes = Lar.Struct([ cube,
        Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
        Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube ])

    V,FV,EV = Lar.struct2lar(threecubes)


        @testset "spaceindex Tests" begin
            @test typeof((V,FV)) == Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
    		@test Lar.spaceindex([.5,.5,.5])((V,FV)) == [5,6,11,12,9]
			@test typeof(Lar.spaceindex([.5,.5,.5])((V,FV))) == Array{Int64,1}
            @test  Lar.spaceindex([0,0,0.])((V,FV)) == [1,19,3,5]
            @test  Lar.spaceindex([0,0,1.])((V,FV)) == [1,19,3,5]
        end

        @testset "rayintersection Tests" begin
            V,(VV,EV,FV,CV) = Lar.simplex(3,true);
    		@test Lar.rayintersection([.333,.333,0])(V,FV,4) == [0.333,0.333,0.3340000000000001]
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "planemap Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "settestpoints Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "testinternalpoint Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "getinternalpoint Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "chainbasis2solids Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "internalpoints Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

        @testset "bool3d Tests" begin
    		@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
			@test bbbb == cccc
        end

end
