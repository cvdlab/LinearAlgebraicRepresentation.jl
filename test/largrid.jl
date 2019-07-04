using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

@testset "LarGrid-0-1 Tests" begin
   @testset "grid_0" begin
      @test Lar.grid_0(1) == [0  1]
      @test Lar.grid_0(2) == [0  1  2]
      @test Lar.grid_0(3) == [0  1  2  3]
   end
   
   @testset "Grid_1 Tests" begin
      @test repr(Lar.grid_1(1)) == "[0; 1]"
      @test Lar.grid_1(2) == [0 1; 1 2]
      @test Lar.grid_1(3) == [0 1 2; 1 2 3]
   end
   
   @testset "Grid Tests" begin
      @test repr(Lar.larGrid(1)(0)) == "[0 1]"
      @test repr(Lar.larGrid(1)(1)) == "[0; 1]"
      @test Lar.larGrid(2)(0) == [0  1  2]
      @test Lar.larGrid(2)(1) == [0 1; 1 2]
      @test Lar.larGrid(3)(0) == [0  1  2  3]
      @test Lar.larGrid(3)(1) == [0  1  2; 1  2  3]
   end
end

@testset "LarVertProd" begin
      vertexDomain(n) = hcat([k for k in 0:n-1]...)
      
      @testset "LarVertProd 1D" begin
            shape = [3]
	      vertLists = [vertexDomain(k+1) for k in shape]
            @test typeof(Lar.larVertProd(vertLists))==Array{Int64,2}
            @test size(Lar.larVertProd(vertLists))==(1, 4)
            @test Lar.larVertProd(vertLists)[:,1]==[0]
            @test Lar.larVertProd(vertLists)[:,4]==[3]
      end

      @testset "LarVertProd 2D" begin
            shape = [3,2]
	      vertLists = [vertexDomain(k+1) for k in shape]
            @test typeof(Lar.larVertProd(vertLists))==Array{Int64,2}
            @test size(Lar.larVertProd(vertLists))==(2, 12)
            @test Lar.larVertProd(vertLists)[:,1]==[0;0]
            @test Lar.larVertProd(vertLists)[:,12]==[3;2]
      end

      @testset "LarVertProd 3D" begin
            shape = [3,2,1]
	      vertLists = [vertexDomain(k+1) for k in shape]
            @test typeof(Lar.larVertProd(vertLists))==Array{Int64,2}
            @test size(Lar.larVertProd(vertLists))==(3, 24)
            @test Lar.larVertProd(vertLists)[:,1]==[0;0;0]
            @test Lar.larVertProd(vertLists)[:,24]==[3;2;1]
      end
end

@testset "Index2addr Tests" begin
   @testset "Shape 1D Tests" begin
      @test Lar.index2addr([10])([0])==1
      @test Lar.index2addr([10])([9])==10
      @test [Lar.index2addr([10])([index]) for index in collect(0:9)]==collect(1:10)
   end

   @testset "Shape 2d Tests" begin
      aa = Lar.cart([[0;1;2],[0;1]])
      bb = Lar.cart([collect(0:9),collect(0:1)])
      cc = Lar.cart([collect(0:2),collect(0:2)])
      dd = "Tuple{Int64,Int64}[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]"
      ee = "Tuple{Int64,Int64}[(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]"
      ff = repr( Lar.cart([collect(0:2),collect(0:1)]) )
      @test [ Lar.index2addr([3,2])(collect(index)) for index in aa ]==collect(1:6)
      @test [ Lar.index2addr([10,2])(collect(index)) for index in bb ] == collect(1:20)
      @test [ Lar.index2addr([3,3])(collect(index)) for index in cc] == collect(1:9)
   end

   @testset "Shape 3d Tests" begin
      aaa = Lar.cart([collect(0:3),collect(0:2),collect(0:1)])
      bbb = Lar.cart([[0;1;2],[0;1],[0;1]])
      ccc = Lar.cart([[0;1],[0;1],[0;1;2]])
      @test Lar.index2addr([3,2,1])([0,0,0]) == 1
      @test [Lar.index2addr([4,3,2])(collect(index)) 
      	for index in aaa ] == collect(1:24)
      @test [Lar.index2addr([3,2,2])(collect(index)) 
      	for index in bbb ] == collect(1:12)
      @test [Lar.index2addr([2,2,3])(collect(index))
      	for index in ccc ] == collect(1:12)
   end
end

@testset "LarCellProd Tests" begin
	 @testset "3D shape" begin
		@test typeof(Lar.larGridSkeleton([3,2,1])(3)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3,2,1])(2)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3,2,1])(1)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3,2,1])(0)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([1,1,1])(3)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([1,1,1])(2)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([1,1,1])(1)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([1,1,1])(0)) == 
		   		Array{Array{Int64,1},1}
	 end
	 @testset "2D shape" begin
		@test typeof(Lar.larGridSkeleton([3,2])(2)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3,2])(1)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3,2])(0)) == 
		   		Array{Array{Int64,1},1}
	 end
	 @testset "1D shape" begin
		@test typeof(Lar.larGridSkeleton([3])(1)) == 
		   		Array{Array{Int64,1},1}
		@test typeof(Lar.larGridSkeleton([3])(0)) == 
		   		Array{Array{Int64,1},1}
	 end
end

@testset "FilterByOrder Tests" begin
      term = "000"
      bit = '0'
      theTerm = Vector{Char}(term)
      @test typeof(theTerm) == Array{Char,1}
      @test parse(Int8,bit) == 0
      @test [parse(Int8,bit) for bit in theTerm] == zeros(3)
      out = hcat([[parse(Int8,bit) for bit in term] for term in 
      	Lar.binaryRange(3)]...)
      @test typeof(out) == Array{Int8,2}
      @test size(out) == (3,8)
      @test repr(out) == "Int8[0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1]"

      @testset "$n" for n in 1:4
            data = [Lar.filterByOrder(n)[k] for (k,el) in enumerate(Lar.filterByOrder(n))]
            @test sum(map(length,data)) == 2^n
      end

      @testset "$n,$k" for n in 1:4, k in 0:n
            @test length(Lar.filterByOrder(n)[k+1]) == binomial(n,k)
      end 
end

@testset "LarImageVerts Tests" "$shape" for shape in [[3,2,1],[3,2],[10,10,10]]
          @test size(Lar.larImageVerts(shape)) == (length(shape),prod(shape .+ 1))
       end

@testset "LarCuboids Tests" begin
      @testset "Simple shape" begin
            shape = [3,2,1]
            cubes = Lar.cuboidGrid(shape,true)
            verts, cells = cubes
            @test typeof(verts) == Array{Float64,2}
            VV,EV,FV,CV = cells
            @testset "$basis" for basis in [VV,EV,FV,CV]
              @test typeof(basis) == Array{Array{Int64,1},1}
            end
      end

      @testset "$shape" for shape in [[3,2,1],[1,1,1],[3,3,10]]
            endShape = shape .+ 1
            cubes = Lar.cuboidGrid(shape,true);
            verts, cells = cubes;
            VV,EV,FV,CV = cells;
            @test size(Lar.cuboidGrid(shape)[1],2) == prod(collect(endShape))
            @test length(Lar.cuboidGrid(shape)[2]) == prod(shape)
            @test typeof(verts) == Array{Float64,2}
	      @test typeof(VV) == Array{Array{Int64,1},1}
	      @test typeof(EV) == Array{Array{Int64,1},1}
	      @test typeof(FV) == Array{Array{Int64,1},1}
	      @test typeof(FV) == Array{Array{Int64,1},1}
      end

      @testset "$shape" for shape in [[3,2],[1,1],[3,10]]
            endShape = shape .+ 1
            @test size(Lar.cuboidGrid(shape)[1],2) == prod(collect(endShape))
            @test length(Lar.cuboidGrid(shape)[2]) == prod(shape)
            cubes = Lar.cuboidGrid(shape,true)
            verts, cells = cubes
            @test typeof(verts) == Array{Float64,2}
            VV,EV,FV = cells
	      @test typeof(VV) == Array{Array{Int64,1},1}
	      @test typeof(EV) == Array{Array{Int64,1},1}
	      @test typeof(FV) == Array{Array{Int64,1},1}
      end

      @testset "$shape" for shape in [[3],[1],[10]]
            endShape = shape .+ 1
            @test size(Lar.cuboidGrid(shape)[1],2) == prod(collect(endShape))
            @test length(Lar.cuboidGrid(shape)[2]) == prod(shape)
            cubes = Lar.cuboidGrid(shape,true)
            verts, cells = cubes
            @test typeof(verts) == Array{Float64,2}
            VV,EV = cells
	      @test typeof(VV) == Array{Array{Int64,1},1}
	      @test typeof(EV) == Array{Array{Int64,1},1}
      end
end

@testset "LarModelProduct Tests" begin
      geom_0 = hcat([[x] for x=0.:10]...)
      topol_0 = [[i,i+1] for i=1:9]
      geom_1 = hcat([[0.],[1.],[2.]]...)
      topol_1 = [[1,2],[2,3]]
      model_0 = (geom_0, topol_0)
      model_1 = (geom_1, topol_1)
      model_2 = Lar.larModelProduct(model_0, model_1)
      model_3 = Lar.larModelProduct(model_2, model_1)

      @testset "Result size" begin
            @test size(model_2[1],2) == size(model_0[1],2) * size(model_1[1],2)
            @test length(model_2[2]) == length(model_0[2]) * length(model_1[2])
            @test size(model_3[1],2) == size(model_2[1],2) * size(model_1[1],2)
            @test length(model_3[2]) == length(model_2[2]) * length(model_1[2])
      end

      @testset "Model 0 type" begin
            m0_V,m0_EV = model_0
            @test typeof(m0_V) == Array{Float64,2}
            @test typeof(m0_EV) == Array{Array{Int64,1},1}
      end
      
      @testset "Model 1 type" begin
            m1_V,m1_EV = model_1
            @test typeof(m1_V) == Array{Float64,2}
            @test typeof(m1_EV) == Array{Array{Int64,1},1}
      end
      
      @testset "Model 2 type" begin
            m2_V,m2_EV = model_2
            @test typeof(m2_V) == Array{Float64,2}
            @test typeof(m2_EV) == Array{Array{Int64,1},1}
      end
      
      @testset "Model 3 type" begin
            m3_V,m3_EV = model_3
            @test typeof(m3_V) == Array{Float64,2}
            @test typeof(m3_EV) == Array{Array{Int64,1},1}
      end

      @testset "Simple product" begin
            modelOne = ([0 0 1 1; 0 1 0 1], Array{Int64,1}[[1, 2, 3, 4]])
            modelTwo = ([0 1], Array{Int64,1}[[1, 2]])
            @test Lar.larModelProduct(modelOne, modelTwo) == Lar.larModelProduct(modelTwo, modelOne)
            @test Lar.larModelProduct(modelTwo, modelTwo)[1] == [0  0  1  1; 0  1  0  1]
            @test Lar.larModelProduct(modelOne,modelOne)[2][1] == Any[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
      end
end
