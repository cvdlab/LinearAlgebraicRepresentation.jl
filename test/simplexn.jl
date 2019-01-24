using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Test

@testset "simplexn" begin
	@testset "extrudeSimplicial" begin
		V = [[0,0] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];
		FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];
		model = (V,FV);
		pattern = repeat([1,2,-3],outer=4);		
		W,FW = Lar.extrudeSimplicial(model, pattern);
		@test typeof(Lar.extrudeSimplicial(model, pattern))==
			Tuple{Array{Int64,2},Array{Array{Int64,1},1}}
		@test typeof(W)==Array{Int64,2}
		@test size(W)==(3,117)
		@test typeof(FW)==Array{Array{Int64,1},1}
		@test length(FW)==144
	end

	@testset "simplexGrid" begin
		#@test Lar.simplexGrid([0])==([0.0], Array{Int64,1}[])
		@test Lar.simplexGrid([1])==([0.0 1.0], Array{Int64,1}[[1, 2]])
		@test Lar.simplexGrid([1,1])==([0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0], 
			Array{Int64,1}[[1, 2, 3], [2, 3, 4]])
		@test typeof(Lar.simplexGrid([1,1,1])[2])==Array{Array{Int64,1},1}
		@test length(Lar.simplexGrid([1,1,1])[2])==6
		V,CV = Lar.simplexGrid([10,10,1])
		@test typeof(V)==Array{Float64,2}
		@test typeof(CV)==Array{Array{Int64,1},1}
		@test size(V)==(3,242)
		@test length(CV)==600
	end

	@testset "simplexFacets" begin
		@test Lar.simplexGrid([1,1])[1]==[0 1 0 1; 0 0 1 1]
		@test Lar.simplexGrid([1,1])[2]==Array{Int64,1}[[1, 2, 3], [2, 3, 4]]
		#=
		@test typeof(Lar.extrudeSimplicial((V,FV), [1])[1])==Array{Float64,2}
		@test typeof(Lar.extrudeSimplicial((V,FV), [1])[2])==Lar.Cells
		@test size(Lar.extrudeSimplicial((V,FV), [1])[1])==(4, 484)
		@test length(Lar.extrudeSimplicial((V,FV), [1])[2])==18
		W,CW = Lar.extrudeSimplicial((V,FV), [1])
		@test typeof(Lar.simplexFacets(CW))==Array{Array{Int64,1},1}
		@test size(Lar.simplexFacets(CW))==(56,)
		@test typeof(Lar.simplexFacets(CW))==Array{Array{Int64,1},1}
		@test size(Lar.simplexFacets(CW))==(56,)
		=#
	end
end


