
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

   @testset "boundary2 Tests" begin
      V,cells = LARLIB.larCuboids([2,1,1],true) 
      VV,EV,FV,CV = cells
      @test typeof(LARLIB.boundary2(FV,EV))==SparseMatrixCSC{Int8,Int64}
      @test size(LARLIB.boundary2(FV,EV))==(length(FV),length(EV)) 
      @test sum(findnz(LARLIB.boundary2(FV,EV))[3])==0
   end

   @testset "chaincomplex Tests" begin
      V,cells = LARLIB.larCuboids([2,1,1],true) 
      VV,EV,FV,CV = cells
      bases,coboundaries = LARLIB.chaincomplex(V,FV,EV)[2:3]
      @test typeof(bases)==Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1},
         Array{Array{Int64,1},1}}
      @test typeof(coboundaries)==Tuple{SparseMatrixCSC{Int8,Int64},
         SparseMatrixCSC{Int8,Int64},SparseMatrixCSC{Int8,Int64}}
      @testset "$basis" for basis in bases
         @test typeof(basis)==Array{Array{Int64,1},1}
      end
      @testset "$coboundary" for coboundary in coboundaries
         @test typeof(coboundary)==SparseMatrixCSC{Int8,Int64}
      end
   end
   
   @testset "collection2model Tests" begin
      @test typeof(collection)==Array{Array{Array,1},1}
   end

   @testset "facetriangulation Tests" begin
      V,(VV,EV,FV,CV) = LARLIB.larCuboids([2,2,1],true)
      W,FW,EW = copy(V),copy(FV),copy(EV)
      collection = Array{Array,1}[]
      for k=1:2
       W,FW,EW = copy(W)+.5,copy(FV),copy(EV)
       append!(collection, [[W,FV,EV]])
      end
      V,FV,EV = LARLIB.collection2model(collection)
      V,bases,coboundaries = LARLIB.chaincomplex(V,FV,EV)
      EV,FV,CV = bases
      cscEV,cscFE,cscCF = coboundaries
      TV = triangulate((1:length(FV),ones(length(FV))),V,FV,cscFE,cscCF)
      @test typeof(TV)==Array{Array{Int64,1},1}
      @test typeof((V,FV))==Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
      @test typeof((1:length(FV),ones(length(FV))))==Tuple{UnitRange{Int64},Array{Float64,1}}
      @testset "$triangle" for triangle in TV
         @test length(triangle)==3
      end
   end
end
