using LARLIB
using Test

@testset "checkStruct" begin
	square=[([[0.575;-0.175] [0.575;0.175] [0.925;-0.175] [0.925;0.175]],
	[[0,1,2,3]])]
	@test LARLIB.checkStruct(square)==length(square[1][1][:,1])
	@test LARLIB.checkStruct(square)==2
end

@testset "apply Tests" begin
	square = ([[0; 0] [0; 1] [1; 0] [1; 1]], [[1, 2, 3,
			4]], [[1,2], [1,3], [2,4], [3,4]])
	@testset "2D" begin
		@testset "apply Translation 2D" begin
			@test typeof(LARLIB.apply(LARLIB.t(-0.5,-0.5),square))== 
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},
			Array{Array{Int64,1},1}}
			@test LARLIB.apply(LARLIB.t(-0.5,-0.5),square)==
			([-0.5 -0.5 0.5 0.5; -0.5 0.5 -0.5 0.5], 
			Array{Int64,1}[[1, 2, 3, 4]], Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4]])
		end
		@testset "apply Scaling 2D" begin
			@test typeof(LARLIB.apply(LARLIB.s(-0.5,-0.5),square))==
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
			@test LARLIB.apply(LARLIB.s(-0.5,-0.5),square)==
			([0.0 0.0 -0.5 -0.5; 0.0 -0.5 0.0 -0.5], Array{Int64,1}[[1, 2, 3, 4]], 
			Array{Int64,1}[[1, 2], [1, 3], [2, 4], [3, 4]])
		end
		@testset "apply Rotation 2D" begin
			@test typeof(LARLIB.apply(LARLIB.r(0),square))==
			Tuple{Array{Float64,2},Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
			@test LARLIB.apply(LARLIB.r(0),square)==square
		end
	end
end

@testset "3D" begin
	cube=LARLIB.cuboid([1,1,1])
	@testset "apply Translation 3D" begin
		@test typeof(LARLIB.apply(LARLIB.t(-0.5,-0.5,-0.5),cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test LARLIB.apply(LARLIB.t(-0.5, -0.5, -0.5),cube) == 
		([-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5; 
		-0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5; -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5], 
		Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]])
	end 
	@testset "apply Scaling 3D" begin
		@test typeof(LARLIB.apply(LARLIB.s(-0.5,-0.5,-0.5),cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test LARLIB.apply(LARLIB.s(-0.5, -0.5, -0.5),cube) == 
		([0.0 0.0 0.0 0.0 -0.5 -0.5 -0.5 -0.5; 
		0.0 0.0 -0.5 -0.5 0.0 0.0 -0.5 -0.5; 0.0 -0.5 0.0 -0.5 0.0 -0.5 0.0 -0.5], 
		Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]])
	end
	@testset "apply Rotation 3D" begin
		@test typeof(LARLIB.apply(LARLIB.r(pi,0,0),cube))==
		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
		@test LARLIB.apply(LARLIB.r(0, 0, 0),cube) ==
		 ([0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0; 
		 0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0], Array{Int64,1}[[1, 2, 3, 4, 5, 6, 7, 8]]) 
	end
end

@testset "cuboid Tests" begin
	square=LARLIB.cuboid([1,1])
	cube=LARLIB.cuboid([1,1,1])
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
	square=LARLIB.cuboid([1,1])
	dim = 2
	structure=LARLIB.Struct([square])
	@test typeof(structure.body) ==
	Array{Tuple{Array{Float64,2},Array{Array{Int64,1},1}},1}
	@test length(LARLIB.traversal(eye(dim+1),[],structure,[]))==length(structure.body)
	@test typeof(LARLIB.traversal(eye(dim+1),[],structure,[]))==Array{Any,1}
end

@testset "Struct Tests" begin
	square=LARLIB.cuboid([1,1])
	@test LARLIB.Struct([square]).body==[square]
	@test LARLIB.Struct([square]).dim==size(square[1],1)
	@test LARLIB.Struct([square]).body[1]==square
end

# to test
#@testset "embedTraversal Tests" begin
#	square = LARLIB.cuboid([1,1])
#	x = LARLIB.Struct([square])
#	@test length(LARLIB.embedTraversal(LARLIB.Struct(),x,1,"New").body[2][1][1])==
#	length(x.body[1][1][1])+1 
#	#in this case n=1, but generally:
#	# length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
#	@test length(embedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])?==
#	length(x.body[1][1][1])+3
#	@test typeof(embedTraversal(deepcopy(x),deepcopy(x),1,"New"))==Struct		
#end

#@testset "embedStruct Tests" begin
#	square = LARLIB.cuboid([1,1])
#	x = LARLIB.Struct([square])	
#	@test length(LARLIB.embedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
#	#in this case n = 1, but generally: 
#	#length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
#	@test length(LARLIB.embedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
#	@test typeof(LARLIB.embedStruct(1)(x))==Struct		
#end

@testset "removeDups Tests" begin
	CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],
		 [4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
	CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]
	@testset "removeDups 3D" begin
		@test length(LARLIB.removeDups(CW1))<= length(CW1)
		@test typeof(LARLIB.removeDups(CW1))==Array{Array{Int64,1},1}
	end
	@testset "removeDups 2D" begin
		@test length(LARLIB.removeDups(CW2))<= length(CW2)
		@test typeof(LARLIB.removeDups(CW2))==Array{Array{Int64,1},1}
	end
end

@testset "struct2lar" begin
   @testset "struct2lar 2D" begin
      square = LARLIB.cuboid([1,1])
      table = LARLIB.apply(LARLIB.t(-0.5,-0.5),square)
      structure = LARLIB.Struct([repeat([table,LARLIB.r(pi/2)],outer=2)...])
      @test typeof(LARLIB.struct2lar(structure))==
      		Tuple{Array{Float64,2},Array{Array{Int64,1},1}}      
	  @test size(LARLIB.struct2lar(structure)[1])==(2, 4)
   end
   @testset "struct2lar 3D" begin
      BV = [[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]]
      V = hcat([[0.0,0.0,0.0],[0.0,0.0,1.0],[0.0,1.0,0.0],[0.0,1.0,1.0],
      [1.0,0.0,0.0],[1.0,0.0,1.0],[1.0,1.0,0.0],[1.0,1.0,1.0]]...)
      block = (V,BV)
      structure = LARLIB.Struct(repeat([block,LARLIB.t(1,0,0)],outer=2));
      @test typeof(LARLIB.struct2lar(structure))==
      	Tuple{Array{Float64,2},Array{Array{Int64,1},1}}
      @test size(LARLIB.struct2lar(structure)[1],1)==3
    end
end
