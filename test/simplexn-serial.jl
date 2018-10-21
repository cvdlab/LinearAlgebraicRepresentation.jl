#####################################################
### From Julia 0.6 to Julia 1.0
### --- SIMPLEXN - serial unit tests ---
### Authors: Fabio Fatelli, Emanuele Loprevite
#####################################################

include("../src/simplexn-serial.jl")

using Test

@testset "Unit Tests" begin
	@testset "larExtrude1" begin
		VOID = [Int64[]],[[0]] # the empty simplicial model
		pattern1 = [1,-1]
		model1 = larExtrude1(VOID,pattern1) # 1D
		#println(model1)
		@test typeof(model1) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model1[1]) == length(pattern1)+1 # num of vertices, no rep
		@test length(model1[2]) == length(filter(x->x>0,pattern1)) # num of 1D-simplices
		
		pattern2 = [1,1,1,-1]
		model2 = larExtrude1(VOID,pattern2) # 1D
		#println(model2)
		@test typeof(model2) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model2[1]) == length(pattern2)+1 # num of vertices, no rep
		@test length(model2[2]) == length(filter(x->x>0,pattern2)) # num of 1D-simplices
		
		pattern3 = [1,1,1,1,1,-1]
		model3 = larExtrude1(VOID,pattern3) # 1D
		#println(model3)
		@test typeof(model3) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model3[1]) == length(pattern3)+1 # num of vertices, no rep
		@test length(model3[2]) == length(filter(x->x>0,pattern3)) # num of 1D-simplices
		
		model1 = larExtrude1(model1,pattern1) # 2D
		#println(model1)
		@test typeof(model1) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model1[1]) == (length(pattern1)+1)^2 # num of vertices, no rep
		@test length(model1[2]) == 2*length(filter(x->x>0,pattern1))^2 # num of 2D-simplices
		
		model2 = larExtrude1(model2,pattern2) # 2D
		#println(model2)
		@test typeof(model2) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model2[1]) == (length(pattern2)+1)^2 # num of vertices, no rep
		@test length(model2[2]) == 2*length(filter(x->x>0,pattern2))^2 # num of 2D-simplices
		
		model3 = larExtrude1(model3,pattern3) # 2D
		#println(model3)
		@test typeof(model3) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model3[1]) == (length(pattern3)+1)^2 # num of vertices, no rep
		@test length(model3[2]) == 2*length(filter(x->x>0,pattern3))^2 # num of 2D-simplices
		
		model1 = larExtrude1(model1,pattern1) # 3D
		#println(model1)
		@test typeof(model1) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model1[1]) == (length(pattern1)+1)^3 # num of vertices, no rep
		@test length(model1[2]) == 6*length(filter(x->x>0,pattern1))^3 # num of 3D-simplices
		
		model2 = larExtrude1(model2,pattern2) # 3D
		#println(model2)
		@test typeof(model2) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model2[1]) == (length(pattern2)+1)^3 # num of vertices, no rep
		@test length(model2[2]) == 6*length(filter(x->x>0,pattern2))^3 # num of 3D-simplices
		
		model3 = larExtrude1(model3,pattern3) # 3D
		#println(model3)
		@test typeof(model3) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model3[1]) == (length(pattern3)+1)^3 # num of vertices, no rep
		@test length(model3[2]) == 6*length(filter(x->x>0,pattern3))^3 # num of 3D-simplices
		
		model = ([[0.0,0],[0,1]],[[1,0],[1,1]]) # 2D with float
		@test typeof(model[1]) == Array{Array{Float64,1},1}
		model = larExtrude1(model,pattern1)
		#println(model)
		@test typeof(model) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
	end
	@testset "larSimplexGrid1" begin
		shape = [3] # 1D
		model1 = larSimplexGrid1(shape)
		#println(model1)
		@test typeof(model1) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test length(model1[1]) == shape[1]+1 # num of vertices, no rep
		@test length(model1[2]) == shape[1] # num of 1D-simplices
		
		shape = [3,5] # 2D
		model2 = larSimplexGrid1(shape)
		#println(model2)
		@test typeof(model2) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		# (shape[1]+1)*(shape[2]+1)
		@test length(model2[1]) == prod(shape)+sum(shape)+1 # num of vertices, no rep
		@test length(model2[2]) == 2*prod(shape) # num of 2D-simplices
		
		shape = [3,5,7] # 3D
		model3 = larSimplexGrid1(shape)
		#println(model3)
		@test typeof(model3) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		# (shape[1]+1)*(shape[2]+1)*(shape[3]+1); num of vertices, no rep
		@test length(model3[1]) == prod(shape)+sum(prod.(collect(combinations(shape,2))))+sum(shape)+1
		@test length(model3[2]) == 6*prod(shape) # num of 3D-simplices
	end
	@testset "larSimplexFacets" begin
		s = larSimplexGrid1([3])[2] # 1D
		#println(s)
		sOut = larSimplexFacets(s)
		#println(sOut)
		@test typeof(sOut) == Array{Array{Int64,1},1}
		d1 = sum([length(s[k]) for k in 1:length(s)])/length(s)
		d2 = sum([length(sOut[k]) for k in 1:length(sOut)])/length(sOut)
		@test d1-1 == d2 # dimension
		# length(s)*binomial(length(s[1]),(length(s[1])-1))
		@test length(sOut) <= length(s)*length(s[1]) # "<=" because no rep
		
		s = larSimplexGrid1([3,5])[2] # 2D
		#println(s)
		sOut = larSimplexFacets(s)
		#println(sOut)
		@test typeof(sOut) == Array{Array{Int64,1},1}
		d1 = sum([length(s[k]) for k in 1:length(s)])/length(s)
		d2 = sum([length(sOut[k]) for k in 1:length(sOut)])/length(sOut)
		@test d1-1 == d2 # dimension
		# length(s)*binomial(length(s[1]),(length(s[1])-1))
		@test length(sOut) <= length(s)*length(s[1]) # "<=" because no rep
		
		s = larSimplexGrid1([3,5,7])[2] # 3D
		#println(s)
		sOut = larSimplexFacets(s)
		#println(sOut)
		@test typeof(sOut) == Array{Array{Int64,1},1}
		d1 = sum([length(s[k]) for k in 1:length(s)])/length(s)
		d2 = sum([length(sOut[k]) for k in 1:length(sOut)])/length(sOut)
		@test d1-1 == d2 # dimension
		# length(s)*binomial(length(s[1]),(length(s[1])-1))
		@test length(sOut) <= length(s)*length(s[1]) # "<=" because no rep
	end
	@testset "quads2tria" begin
		modIn = ([[0,0,0],[0,1,0],[1,0,0],[1,1,0]],[[0,1,2,3]]) # 2D
		modOut = quads2tria(modIn)
		#println(modOut)
		@test typeof(modIn) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test typeof(modOut) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		edges = 2*sum([length(modIn[2][k]) for k in 1:length(modIn[2])])
		EulerChar = length(modOut[1])-edges+length(modOut[2]) # V - E + F
		@test EulerChar == 1
		
		modIn = ([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0.5,1.5,0]],[[0,1,2,3,4]]) # 2D
		modOut = quads2tria(modIn)
		#println(modOut)
		@test typeof(modIn) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		@test typeof(modOut) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		edges = 2*sum([length(modIn[2][k]) for k in 1:length(modIn[2])])
		EulerChar = length(modOut[1])-edges+length(modOut[2]) # V - E + F
		@test EulerChar == 1
		
		modIn = ([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,1]],[[0,1,2,3,4,5,6,7]]) # 3D
		modOut = quads2tria(modIn)
		#println(modOut)
		@test typeof(modIn) == Tuple{Array{Array{Int64,1},1},Array{Array{Int64,1},1}}
		@test typeof(modOut) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		edges = 2*sum([length(modIn[2][k]) for k in 1:length(modIn[2])])
		EulerChar = length(modOut[1])-edges+length(modOut[2]) # V - E + F
		@test EulerChar == 1
		
		modIn = ([[0,0,0],[0,1,0],[1,0,0],[1,1,0],[0,0,1],[0,1,1],[1,0,1],[1,1,1],[0.5,0.5,1.5]],[[0,1,2,3,4,5,6,7,8]]) # 3D
		modOut = quads2tria(modIn)
		#println(modOut)
		@test typeof(modIn) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		@test typeof(modOut) == Tuple{Array{Array{Float64,1},1},Array{Array{Int64,1},1}}
		edges = 2*sum([length(modIn[2][k]) for k in 1:length(modIn[2])])
		EulerChar = length(modOut[1])-edges+length(modOut[2]) # V - E + F
		@test EulerChar == 1
	end
end
