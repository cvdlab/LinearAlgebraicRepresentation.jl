using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using DataStructures,IntervalTrees


@testset "2D containment tests" begin

(V, EV) = ([0.43145 0.596771 0.758062 1.0 0.778226 0.919353 0.879033 0.806447 0.778226 0.709677 0.596771 0.262094 0.322578 0.0 0.2379 0.161291 0.467739 0.429435 0.627999 0.627999 0.383062 0.694833 0.653221 0.544027 0.778226 0.848789 0.750707 0.627999 0.694833 0.806447; -0.0163938 0.22521 0.104412 0.325182 0.629266 0.683418 0.820882 0.725074 0.845873 0.75215 1.0 0.820882 0.629266 0.385151 0.43765 0.246033 0.466811 0.629266 0.704244 0.507207 0.275195 0.683418 0.43765 0.291323 0.199264 0.43765 0.497413 0.341841 0.259902 0.364484], Array{Int64,1}[[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 12], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18], [18, 19], [19, 20], [20, 21], [21, 1], [22, 23], [24, 23], [24, 25], [25, 26], [26, 22], [27, 28], [28, 29], [29, 30], [30, 27]])

classify = Lar.pointInPolygonClassification(V,EV)
queryPoint = [0.5,0.5]

	@testset "crossingTest Tests" begin
		@test Lar.crossingTest(0, 0, 0., 0)::Number == 0.5
		@test Lar.crossingTest(0, 0, 0.5, 0)::Number == 1.0
		@test Lar.crossingTest(0, 0, 0.5, 0)::Number == 1.0
		@test Lar.crossingTest(1, 0, 0.5, 0) == 1.0
		@test Lar.crossingTest(1, 1, 0.5, 0) == 1.0
		@test Lar.crossingTest(1, 1, 0.5, 1) == 0
	end
	
	@testset "setTile Tests" begin
		x,y = 0.5,0.75
		xmin,xmax,ymin,ymax = x,x,y,y
		box = [ymax,ymin,xmax,xmin]
		tilecode = Lar.setTile(box)
	
		@test Lar.setTile isa Function
		@test typeof(box)==Array{Float64,1}
		@test Lar.setTile(box) isa Function
		@test tilecode([.5,.5])==2
		@test tilecode([-.5,.5])==10
		@test tilecode([.5,-.5])==2
		@test tilecode([-.5,-.5])==10
		@test tilecode([.5,.95])==1
	end
	
	@testset "pointInPolygonClassification Tests" begin
		@test Lar.pointInPolygonClassification(V,EV) isa Function
		@test pnt = [0.5,0.5] isa Array{Float64,1}
		@test classify(queryPoint)=="p_out"
		@test classify([0.5,0.75])=="p_in"
		@test classify([1.5,0.75])=="p_out"
		@test typeof(classify(queryPoint))==String
	end
end



@testset "Biconnected components" begin

(V, (VV, EV, FV)) = Lar.cuboidGrid([3, 3], true)

	@testset "verts2verts data Tests" begin
		@test Lar.cuboidGrid([3,3],true) == ([0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 3.0 3.0 3.0 3.0; 0.0 1.0 2.0 3.0 0.0 1.0 2.0 3.0 0.0 1.0 2.0 3.0 0.0 1.0 2.0 3.0], ( [[[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16]], [[1, 2],[2,3],[3,4],[5,6],[6,7],[7,8],[9,10],[10,11],[11,12],[13,14],[14, 15],[15,16],[1,5],[2,6],[3,7],[4,8],[5,9],[6,10],[7,11],[8,12],[9,13], [10,14],[11,15],[12,16]], [[1,2,5,6],[2,3,6,7],[3,4,7,8],[5,6,9,10], [6,7,10,11],[7,8,11,12],[9,10,13,14],[10,11,14,15],[11,12,15,16]]] ) )
		@test size(V,2)==length(VV)
		@test length(VV)==16
		@test length(EV)==24
		@test length(FV)==9
		@test length(VV)-length(EV)+length(FV)==1
		@test Lar.verts2verts(EV::Lar.Cells)==[[2,5],[1,3,6],[2,4,7],[3,8],[1,
		6,9],[2,5,7,10],[3,6,8,11],[4,7,12],[5,10,13],[6,9,11,14],[7,10,12,
		15],[8,11,16],[9,14],[10,13,15],[11,14,16],[12,15]]
		@test Lar.verts2verts(FV::Lar.Cells)==[[2,5,6],[1,3,5,6,7],[2,4,6,7,
		8],[3,7,8],[1,2,6,9,10],[1,2,3,5,7,9,10,11],[2,3,4,6,8,10,11,12],[3,
		4,7,11,12],[5,6,10,13,14],[5,6,7,9,11,13,14,15],[6,7,8,10,12,14,15,
		16],[7,8,11,15,16],[9,10,14],[9,10,11,13,15],[10,11,12,14,16],[11,12,15]]
	end
	
	@testset "DFV_visit Tests" begin
		@test true
	end
	
	@testset "outputComp Tests" begin
		@test Lar.outputComp isa Function
		@test true
	end
	
	@testset "biconnectedComponent Tests" begin
		#(V, EV) = ([0.0 0.97721 0.97721 0.724048 0.724048 0.258225 0.258225 0.660757 0.660757 0.0; 1.0 1.0 0.0 0.0 0.934178 0.934178 0.346836 0.346836 0.0 0.0], [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 1]])
		#V,EVs = Lar.biconnectedComponent((V,EV))
		#@test sort(map(sort,EVs[1]))==sort(map(sort,EV))
		#@test Lar.biconnectedComponent((V,EV)) == ([0.0 0.97721 0.97721 0.724048 0.724048 0.258225 0.258225 0.660757 0.660757 0.0; 1.0 1.0 0.0 0.0 0.934178 0.934178 0.346836 0.346836 0.0 0.0], Any[Array{Int64,1}[[1, 10], [9, 10], [8, 9], [7, 8], [6, 7], [5, 6], [4, 5], [3, 4], [2, 3], [1, 2]]])
		#@test length(EVs)==1
		#@test typeof(EVs[1])==Array{Array{Int64,1},1}
		#@test typeof(EVs)==Array{Any,1}
	end
end


@testset "Refactoring spaceindex tests" begin

	# 2x2x2 cuboidal grid for 1-, 2-, and 3-dim tests
	V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
	W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))

	function test_bboxes(bboxes)
		# initialize accumulator
		accumulator = BitArray{1}()
		for k=1:size(bboxes[1],1)
			push!(accumulator, true)
		end
		# testing data 
		for h=1:length(bboxes)
			accumulator = (bboxes[h][:,1] .< bboxes[h][:,2]) .& accumulator
		end
		return (&)(accumulator...)
	end
	
	@testset "boundingbox Tests" begin
			
		@testset "Edge tests" begin # 
			cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
			bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
		@testset "Face tests" begin # 
			cellpoints = [ W[:,FV[k]]::Lar.Points for k=1:length(FV) ]
			bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
		@testset "Cell tests" begin # 
			cellpoints = [ W[:,CV[k]]::Lar.Points for k=1:length(CV) ]
			bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
			@test true == test_bboxes(bboxes)
		end
	end

	@testset "coordintervals Tests" begin
	# 2x2x2 cuboidal grid for 1-, 2-, and 3-dim tests
	V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
	W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))
			
		cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		@testset "Edge tests" begin # 
			@test typeof(Lar.coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
		cellpoints = [ W[:,FV[k]]::Lar.Points for k=1:length(FV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		@testset "Face tests" begin # 
			@test typeof(Lar.coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
		cellpoints = [ W[:,CV[k]]::Lar.Points for k=1:length(CV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		@testset "Cell tests" begin # 
			@test typeof(Lar.coordintervals(1,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(2,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
			@test typeof(Lar.coordintervals(3,bboxes)) == 
				OrderedDict{Array{Float64,1}, Array{Int64,1}}
		end
	end
	
	@testset "boxcovering Tests" begin
		V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)
		W,_ = Lar.apply(Lar.r(1,1,pi/6),(V,[VV,EV,FV,CV]))
		cellpoints = [ W[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		dict = Lar.coordintervals(1,bboxes)
		@test typeof(dict) == OrderedDict{Array{Float64,1},Array{Int64,1}}
		@test length(Lar.coordintervals(1,bboxes)) == 54
		@test length(Lar.coordintervals(2,bboxes)) == 54
		@test length(Lar.coordintervals(3,bboxes)) == 54

		V,(VV,EV,FV) = Lar.cuboidGrid([2,1],true)
		cellpoints = [ V[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		@test bboxes == [[0.0 0.0; 0.0 1.0],
        [1.0 1.0; 0.0 1.0],
        [2.0 2.0; 0.0 1.0],
        [0.0 1.0; 0.0 0.0],
        [0.0 1.0; 1.0 1.0],
        [1.0 2.0; 0.0 0.0],
        [1.0 2.0; 1.0 1.0]]
        xboxdict = Dict(
         [0.0, 0.0] => [1],
         [1.0, 1.0] => [2],
         [2.0, 2.0] => [3],
         [0.0, 1.0] => [4, 5],
         [1.0, 2.0] => [6, 7])
        @test xboxdict == Lar.coordintervals(1,bboxes)
		xs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in xboxdict
			xs[tuple(key...)] = boxset
		end
       @test typeof(xs) ==
		IntervalTrees.IntervalBTree{Float64,
		IntervalValue{Float64,Array},64}
	end
	
	@testset "Refactoring spaceindex tests" begin
		V,(VV,EV,FV) = Lar.cuboidGrid([2,1],true)
		EV = [[1, 2], [3, 4], [5, 6], [1, 3], [2, 4], [3, 5], [4, 6]]
		cellpoints = [ V[:,EV[k]]::Lar.Points for k=1:length(EV) ]
		bboxes = [hcat(Lar.boundingbox(cell)...) for cell in cellpoints]
		xboxdict = Lar.coordintervals(1,bboxes)
		yboxdict = Lar.coordintervals(2,bboxes)
		xs = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in xboxdict
			xs[tuple(key...)] = boxset
		end
		ys = IntervalTrees.IntervalMap{Float64, Array}()
		for (key,boxset) in yboxdict
			ys[tuple(key...)] = boxset
		end
		xcovers = Lar.boxcovering(bboxes, 1, xs)
		ycovers = Lar.boxcovering(bboxes, 2, ys)
		covers = [intersect(pair...) for pair in zip(xcovers,ycovers)]
		
		@test covers == Array{Int64,1}[[1, 4, 5], [4, 5, 2, 6, 7], [6, 7, 3], 
			[1, 4, 2, 6], [1, 5, 2, 7], [4, 2, 6, 3], [5, 2, 7, 3]]
	end
end


@testset "Refactoring fragmentlines" begin

	@testset "linefragments Tests" begin
		V = hcat([[0.,0],[1,0],[1,1],[0,1],[2,1]]...);
		EV = [[1,2],[2,3],[3,4],[4,1],[1,5]];
		@test Lar.spaceindex((V,EV)) == 
		[[4, 5, 2], [1, 3, 5], [4, 5, 2], [1, 3, 5], [4, 1, 3, 2]]
		Sigma = [[4, 5, 2], [1, 3, 5], [4, 5, 2], [1, 3, 5], [4, 1, 3, 2]]
		@test Lar.linefragments(V,EV,Sigma) ==
		[[0.0, 1.0], [0.0, 0.5, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 0.5, 1.0]]
	end
	
	@testset "intersection Tests" begin
		line1 = [[0.0, 0.0], [1.0, 2.0]]
		line2 = [[2.0, 0.0], [0.0, 3.0]]
 		@test Lar.intersection(line1,line2)==(0.8571428571428571, 0.5714285714285714)
		line1 = [[0.0, 0.0], [-2.0, 0.0]]
		line2 = [[0.0, 0.0], [0.0, 3.0]]
		@test Lar.intersection(line1,line2)==(0.0, 0.0)
		line2 = [[0.0, 3.0],[0.0, 0.0]]
		@test Lar.intersection(line1,line2)==(0.0, 1.0)
		line1 = [[-2.0, 0.0],[0.0, 0.0]]
		@test Lar.intersection(line1,line2)==(1.0, 1.0)
		@test Lar.intersection(line1,line1)==nothing
	end
	
	@testset "fragmentlines Tests" begin
		V = hcat([[0.,0],[1,0],[1,1],[0,1],[2,1]]...);
		EV = [[1,2],[2,3],[3,4],[4,1],[1,5]];
		W,EW = Lar.fragmentlines((V,EV))
		@test W == [0.0  1.0  1.0  1.0  0.0  2.0; 0.0  0.0  0.5  1.0  1.0  1.0]
		@test EW == [[1, 2],[2, 3],[3, 4],[4, 5],[5, 1],[1, 3],[3, 6]]
	end
end


@testset "Refactoring pipeline 2" begin

	@testset "bbbbbbb Tests" begin
		@test true
	end
	
	@testset "bbbbbbb Tests" begin
		@test true
	end
	
	@testset "bbbbbbb Tests" begin
		@test true
	end
end

