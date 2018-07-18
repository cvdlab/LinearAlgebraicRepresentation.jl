
using Base.Test

@testset "checkStruct" begin
	square=[([[0.575;-0.175] [0.575;0.175] [0.925;-0.175] [0.925;0.175]],
	[[0,1,2,3]])]
	@test checkStruct(square)==length(square[1][1][:,1])
	@test checkStruct(square)==2
end


@testset "apply Tests" begin
	square = ([[0; 0] [0; 1] [1; 0] [1; 1]], [[1, 2, 3,
			4]], [[1,2], [1,3], [2,4], [3,4]])
	@testset "2D" begin
		@testset "apply Translation 2D" begin
			@test typeof(apply(t(-0.5,-0.5))(square))== 
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},
			Array{Array{Int64,1},1}}
			@test apply(t(-0.5,-0.5))(square)==([-0.5 -0.5 0.5 0.5; -0.5 0.5 -0.5 0.5], 
			Array{Int64,1}[[1, 2, 3, 4]], Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4]])
		end
		@testset "apply Scaling 2D" begin
			@test typeof(apply(s(-0.5,-0.5))(square))==
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
			@test apply(s(-0.5,-0.5))(square)==
			([0.0 0.0 -0.5 -0.5; 0.0 -0.5 0.0 -0.5], Array{Int64,1}[[1, 2, 3, 4]], 
			Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4]])
		end
		@testset "apply Rotation 2D" begin
			@test typeof(apply(r(0))(square))==
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
			@test apply(r(0))(square)==square
		end
	end
end

@testset "3D" begin
	@testset "apply Translation 3D" begin
		@test typeof(apply(t(-0.5,-0.5,-0.5))(cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test (apply(t(-0.5, -0.5, -0.5)))(cube) == ([-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5; 
		-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5; -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5], 
		Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]])
	end 
	@testset "apply Scaling 3D" begin
		@test typeof(apply(s(-0.5,-0.5,-0.5))(cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test (apply(s(-0.5, -0.5, -0.5)))(cube) == ([0.0 0.0 0.0 0.0 -0.5 -0.5 -0.5 -0.5; 
		0.0 0.0 -0.5 -0.5 0.0 0.0 -0.5 -0.5; 0.0 -0.5 0.0 -0.5 0.0 -0.5 0.0 -0.5], 
		Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]])


	end
	@testset "apply Rotation 3D" begin
		@test typeof(apply(r(pi,0,0))(cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test (apply(r(0, 0, 0)))(cube) ==
		 ([0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0; 
		 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0], Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]]) 
	end
end




@testset "cuboid Tests" begin
	square=cuboid([1,1])
	cube=cuboid([1,1,1])
	@testset "cuboid Tests 2D" begin
		@test typeof(square)==Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test length(square)==2
		@test size(square[1])==(2, 4)
	end
	@testset "cuboid Tests 3D" begin
		@test typeof(cube)==Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test length(cube)==2
		@test size(cube[1])==(3, 8)
	end
end




@testset "traversal" begin
   square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
   @everywhere structure=Struct([square])
   @everywhere dim=checkStruct(structure.body)
   @test length(traversal(eye(dim+1),[],structure,[]))==length(structure.body)
   @test typeof(traversal(eye(dim+1),[],structure,[]))==ArrayAny,1
end


square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

@testset "Struct Tests" begin
	@test Struct([square]).body==[square]
	@test Struct([square]).dim==length(square[1][1])
	@test Struct([square]).cuboid==[[0,0],[1,1]]
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
@testset "pStruct Tests" begin
@test pStruct([square]).body==[square]
	@test pStruct([square]).dim==length(square[1][1])
	@test pStruct([square]).cuboid==[[0,0],[1,1]]
	@test pStruct([square]).category=="feature"
	@test pStruct([square],"quadrato").name=="quadrato"
end



square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])

@testset "embedTraversal Tests" begin
	@test length(embedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1])?==length(x.body[1][1][1])+1 
	#in this case n=1, but generally:
	# length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
	@test length(embedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])?==length(x.body[1][1][1])+3
	@test typeof(embedTraversal(deepcopy(x),deepcopy(x),1,"New"))==Struct		
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])

@testset "pembedTraversal Tests" begin
	@test length(pembedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1])?==length(x.body[1][1][1])+1 
	#in this case n=1, but generally:
	# length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
	@test length(pembedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])==?length(x.body[1][1][1])+3
	@test typeof(pembedTraversal(deepcopy(x),deepcopy(x),1,"New"))==pStruct		
end



square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])	
@testset "embedStruct Tests" begin
	@test length(embedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
	#in this case n = 1, but generally: 
	#length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
	@test length(embedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
	@test typeof(embedStruct(1)(x))==Struct		
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])	
@testset "pembedStruct Tests" begin	
	@test length(pembedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
	#in this case n = 1, but generally: 
	#length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
	@test length(pembedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
	@test typeof(pembedStruct(1)(x))==pStruct		
end


CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],
     [4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "removeDups Tests" begin
   @testset "removeDups 3D" begin
      @test length(removeDups(CW1))<= length(CW1)
      @test typeof(removeDups(CW1))==ArrayArrayInt64,1,1
   end
   @testset "removeDups 2D" begin
      @test length(removeDups(CW2))<= length(CW2)
      @test typeof(removeDups(CW2))==ArrayArrayInt64,1,1
   end
end

CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],
     [4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "premoveDups Tests" begin
   @testset "premoveDups 3D" begin
      @test length(premoveDups(CW1))<= length(CW1)
      @test typeof(premoveDups(CW1))==ArrayArrayInt64,1,1
   end
   @testset "premoveDups 2D" begin
      @test length(premoveDups(CW2))<= length(CW2)
      @test typeof(premoveDups(CW2))==ArrayArrayInt64,1,1
    end
end



@testset "struct2lar" begin
   @testset "struct2lar 2D" begin
      square=([[0, 0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
      table=apply(t(-0.5,-0.5))(square)
      structure=Struct([repeat([table,r(pi/2)],outer=2)...])
      @test typeof(struct2lar(structure))==TupleArrayAny,1,
      ArrayArrayAny,1,1
      @test length(struct2lar(structure)[1][1])==2
   end
   @testset "struct2lar 3D" begin
       BV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
      V=[[0 0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],
         [1,0,1],[1,1,0],[1,1,1]]
      block=[V,BV]
      structure=Struct(repeat([block,t(1,0,0)],outer=2));
      @test typeof(struct2lar(structure))==TupleArrayAny,1,
      ArrayArrayAny,1,1
      @test length(struct2lar(structure)[1][1])==3
    end
end

