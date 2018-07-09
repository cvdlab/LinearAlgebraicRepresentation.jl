@testset "checkStruct" begin{
 list=([[0.575,-0.175],[0.575,0.175],[0.925,-0.175],[0.925,0.175]],{
 [[0,1,2,3]])}
@test checkStruct(list)==length(list[1][1][1])
@test typeof(checkStruct(list))==Int
end


@testset "larApply Tests" begin
   @testset "2D" begin
      @testset "larApply Translation 2D" begin
      @test typeof(larApply(t(-0.5,-0.5))(square))==Tuple{Array{Array{
                 Float64,N} where N,1},Array{Array{Int64,1},1}}
      @test larApply(t(-0.5,-0.5))(square)==([[-0.5,-0.5],[-0.5,0.5],
           [0.5,-0.5],[0.5,0.5]],[[0,1,2,3]])
   end
      @testset "larApply Scaling 2D" begin
         @test typeof(larApply(s(-0.5,-0.5))(square))==Tuple{Array{Array{
                  Float64,N} where N,1},Array{Array{Int64,1},1}}
         @test larApply(s(-0.5,-0.5))(square)==([[0.0,0.0],[0.0,-0.5],
                [-0.5,0.0],[-0.5,-0.5]],[[0,1,2,3]])
      end
      @testset "larApply Rotation 2D" begin
         @test typeof(larApply(r(0))(square))==Tuple{Array{Array{
                 Float64,N}  where N,1},Array{Array{Int64,1},1}}
         @test larApply(r(0))(square)==square
      end
   end

   @testset "3D" begin
      @testset "larApply Translation 3D" begin
         @test typeof(larApply(t(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{
               Array{Float64,N} where N,1},Array{Array{Int64,1},1}}
         @test larApply(t(-0.5,-0.5,-0.5))(cubes)==([[-0.5,-0.5,-0.5],
               [-0.5,-0.5,0.5],[-0.5,0.5,-0.5],[-0.5,0.5,0.5],
               [0.5,-0.5,-0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],
               [0.5,0.5,0.5]],[[0,1,2,3,4,5,6,7]])
      end 
	
      @testset "larApply Scaling 3D" begin
         @test typeof(larApply(s(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{
               Array{Float64,N} where N,1},Array{Array{Int64,1},1}}
         @test larApply(s(-0.5,-0.5,-0.5))(cubes)==([[0.0,0.0,0.0],
                 [0.0,0.0,-0.5],[0.0,-0.5,0.0],[0.0,-0.5,-0.5],
                 [-0.5,0.0,0.0],[-0.5,0.0,-0.5],[-0.5,-0.5,0.0],
                 [-0.5,-0.5,-0.5]],[[0,1,2,3,4,5,6,7]])
      end
	
      @testset "larApply Rotation 3D" begin
         @test typeof(larApply(r(pi,0,0))(cubes))==Tuple{Array{Array{
                Float64,N} where N,1},Array{Array{Int64,1},1}}


@testset "plarApply Tests" begin
   @testset "2D" begin	
      @testset "plarApply Translation 2D" begin
         @test typeof(plarApply(t(-0.5,-0.5))(square))==Tuple{Array{
                    Array{Float64,1},1},Array{Array{Int64,1},1}}
         @test plarApply(t(-0.5,-0.5))(square)==([[-0.5,-0.5],[-0.5,0.5],
                    [0.5,-0.5],[0.5,0.5]],[[0,1,2,3]])
      end 
	
      @testset "plarApply Scaling 2D" begin
         @test typeof(plarApply(s(-0.5,-0.5))(square))==Tuple{Array{
                   Array{Float64,1},1},Array{Array{Int64,1},1}}
         @test plarApply(s(-0.5,-0.5))(square)==([[0.0,0.0],[0.0,-0.5],
                   [-0.5,0.0],[-0.5,0.5]],[[0,1,2,3]])
      end
	
      @testset "plarApply Rotation 2D" begin
         @test typeof(plarApply(r(0))(square))==Tuple{Array{Array{
                   Float64,1},1},Array{Array{Int64,1},1}}
         @test plarApply(r(0))(square)==square
      end
   end
  @testset "3D" begin
    @testset "plarApply Translation 3D" begin
      @test typeof(plarApply(t(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{
                 Array{Float64,1},1},Array{Array{Int64,1},1}}
      @test plarApply(t(-0.5,-0.5,-0.5))(cubes)==([[-0.5,-0.5,-0.5],
                 [-0.5,-0.5,0.5],[-0.5,0.5,-0.5],[-0.5,0.5,0.5],
                 [0.5,-0.5,-0.5],[0.5,-0.5,0.5],[0.5,0.5,-0.5],
                 [0.5,0.5,0.5]],[[0,1,2,3,4,5,6,7]])
    end 
    @testset "plarApply Scaling 3D" begin
      @test typeof(plarApply(s(-0.5,-0.5,-0.5))(cubes))==Tuple{Array{
                 Array{Float64,1},1},Array{Array{Int64,1},1}}
      @test plarApply(s(-0.5,-0.5,-0.5))(cubes)==([[0.0,0.0,0.0],
                 [0.0,0.0,-0.5],[0.0,-0.5,0.0],[0.0,-0.5,-0.5],
                 [-0.5,0.0,0.0],[-0.5,0.0,-0.5],[-0.5,-0.5,0.0],
                 [-0.5,-0.5,-0.5]],[[0,1,2,3,4,5,6,7]])
    end
    @testset "plarApply Rotation 3D" begin
      @test typeof(plarApply(r(pi,0.0,0.0))(cubes))==Tuple{Array{Array{
                    Float64,1},1},Array{Array{Int64,1},1}}


square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],
       [1,1,1]],[[[0],[1],[2],[3],[4],[5],[6],[7]],[[0,1],[2,3],[4,5],
       [6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]],
       [[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],
       [1,3,5,7]],[[0,1,2,3,4,5,6,7]]])

@testset "box Tests" begin
  @testset "box Tests 2D" begin
    @test typeof(box(square))==Array{Array,1}
    @test length(box(square))==2
    @test length(box(square)[1])==2
  end
  @testset "box Tests 3D" begin
    @test typeof(box(cubes))==Array{Array,1}
    @test length(box(cubes))==2
    @test length(box(cubes)[1])==3
  end
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
cubes=([[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],
       [1,1,1]],[[[0],[1],[2],[3],[4],[5],[6],[7]],[[0,1],[2,3],[4,5],
       [6,7],[0,2],[1,3],[4,6],[5,7],[0,4],[1,5],[2,6],[3,7]],
       [[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],
       [1,3,5,7]],[[0,1,2,3,4,5,6,7]]])

@testset "pbox Tests" begin
	@testset "pbox Tests 2D" begin
		@test typeof(pbox(square))==Array{Array,1}
		@test length(pbox(square))==2
		@test length(pbox(square)[1])==2
	end
	@testset "pbox Tests 3D" begin
		@test typeof(pbox(cubes))==Array{Array,1}
		@test length(pbox(cubes))==2
		@test length(pbox(cubes)[1])==3
	end
end



@testset "traversal" begin
   square=([[0, 0], [0, 1], [1, 0], [1, 1]], [[0, 1, 2, 3]])
   @everywhere structure=Struct([square])
   @everywhere dim=checkStruct(structure.body)
   @test length(traversal(eye(dim+1),[],structure,[]))==length(structure.body)
   @test typeof(traversal(eye(dim+1),[],structure,[]))==Array{Any,1}
end


square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])

@testset "Struct Tests" begin{
	@test Struct([square]).body==[square]
	@test Struct([square]).dim==length(square[1][1])
	@test Struct([square]).box==[[0,0],[1,1]]}
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
@testset "pStruct Tests" begin{
@test pStruct([square]).body==[square]{
	@test pStruct([square]).dim==length(square[1][1])
	@test pStruct([square]).box==[[0,0],[1,1]]
	@test pStruct([square]).category=="feature"
	@test pStruct([square],"quadrato").name=="quadrato"}
end



square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])

@testset "embedTraversal Tests" begin{
	@test length(embedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1]){?==length(x.body[1][1][1])+1 }}
	#in this case n=1, but generally:
	# length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
	@test length(embedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1]){?==length(x.body[1][1][1])+3}}
	@test typeof(embedTraversal(deepcopy(x),deepcopy(x),1,"New"))==Struct	}	
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])

@testset "pembedTraversal Tests" begin{
	@test length(pembedTraversal(deepcopy(x),deepcopy(x),1,"New").body[2][1][1]){?==length(x.body[1][1][1])+1 }}
	#in this case n=1, but generally:
	# length(length(embedTraversal(x,x,1,"New")=length(x.body[1][1][1])+n
	@test length(pembedTraversal(deepcopy(x),deepcopy(x),3,"New").body[2][1][1])=={?length(x.body[1][1][1])+3}}
	@test typeof(pembedTraversal(deepcopy(x),deepcopy(x),1,"New"))==pStruct	}	
end



square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=Struct([square])	
@testset "embedStruct Tests" begin{
	@test length(embedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
	#in this case n = 1, but generally: 
	#length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
	@test length(embedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
	@test typeof(embedStruct(1)(x))==Struct	}	
end

square=([[0,0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
x=pStruct([square])	
@testset "pembedStruct Tests" begin{	
	@test length(pembedStruct(1)(x).body[1][1][1])==length(x.body[1][1][1])+1 
	#in this case n = 1, but generally: 
	#length(embedStruct(n)(x).body[1][1][1])=length(x.body[1][1][1])+n
	@test length(pembedStruct(3)(x).body[1][1][1])==length(x.body[1][1][1])+3
	@test typeof(pembedStruct(1)(x))==pStruct	}	
end


CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],
     [4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "removeDups Tests" begin
   @testset "removeDups 3D" begin
      @test length(removeDups(CW1))<= length(CW1)
      @test typeof(removeDups(CW1))==Array{Array{Int64,1},1}
   end
   @testset "removeDups 2D" begin
      @test length(removeDups(CW2))<= length(CW2)
      @test typeof(removeDups(CW2))==Array{Array{Int64,1},1}
   end
end

CW1=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7],
     [4,5,6,7],[8,9,10,11],[4,5,8,9],[6,7,10,11],[4,6,8,10],[5,7,9,11]]
CW2=[[0,1,2,3],[4,5,6,7],[8,9,10,11],[12,13,14,15],[16,17,18,19]]

@testset "premoveDups Tests" begin
   @testset "premoveDups 3D" begin
      @test length(premoveDups(CW1))<= length(CW1)
      @test typeof(premoveDups(CW1))==Array{Array{Int64,1},1}
   end
   @testset "premoveDups 2D" begin
      @test length(premoveDups(CW2))<= length(CW2)
      @test typeof(premoveDups(CW2))==Array{Array{Int64,1},1}
    end
end



@testset "struct2lar" begin
   @testset "struct2lar 2D" begin
      square=([[0, 0],[0,1],[1,0],[1,1]],[[0,1,2,3]])
      table=larApply(t(-0.5,-0.5))(square)
      structure=Struct([repeat([table,r(pi/2)],outer=2)...])
      @test typeof(struct2lar(structure))==Tuple{Array{Any,1},
      Array{Array{Any,1},1}}
      @test length(struct2lar(structure)[1][1])==2
   end
   @testset "struct2lar 3D" begin
       BV=[[0,1,2,3],[4,5,6,7],[0,1,4,5],[2,3,6,7],[0,2,4,6],[1,3,5,7]]
      V=[[0 0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],
         [1,0,1],[1,1,0],[1,1,1]]
      block=[V,BV]
      structure=Struct(repeat([block,t(1,0,0)],outer=2));
      @test typeof(struct2lar(structure))==Tuple{Array{Any,1},
      Array{Array{Any,1},1}}
      @test length(struct2lar(structure)[1][1])==3
    end
end

