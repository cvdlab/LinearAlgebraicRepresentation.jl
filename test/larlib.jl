
@testset "Boundary operators Tests" begin
   @testset "characteristicMatrix" begin
      V,cells = LARLIB.larCuboids([2,1,1],true) 
      VV,EV,FV,CV = cells
      @test typeof(V)==Array{Float64,2}
      @test typeof(cells)==Array{Array{Array{Int64,1},1},1}
      @testset "$rel" for rel in cells
         @test typeof(rel)==Array{Array{Int64,1},1}
         @test typeof(LARLIB.characteristicMatrix(rel))==SparseMatrixCSC{Int8,Int64}
      end
   end

   @testset "boundary1 Tests" begin
      V,cells = LARLIB.larCuboids([2,1,1],true) 
      VV,EV,FV,CV = cells
      @testset "$rel" for rel in cells
         @test typeof(rel)==Array{Array{Int64,1},1}
         @test typeof(LARLIB.boundary1(rel))==SparseMatrixCSC{Int8,Int64}
         @test size(LARLIB.boundary1(EV))==(size(V,2),length(EV)) 
      end
   end

   @testset "uboundary2 Tests" begin
      V,cells = LARLIB.larCuboids([2,1,1],true) 
      VV,EV,FV,CV = cells
      @test typeof(LARLIB.uboundary2(FV,EV))==SparseMatrixCSC{Int8,Int64}
      @test size(LARLIB.uboundary2(FV,EV))==(length(FV),length(EV)) 
   end
end
