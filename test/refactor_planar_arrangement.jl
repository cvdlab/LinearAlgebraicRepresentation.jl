using Test
using SparseArrays
using LinearAlgebra
using LinearAlgebraicRepresentation
using LinearAlgebraicRepresentation.Arrangement
Lar = LinearAlgebraicRepresentation
Lara = LinearAlgebraicRepresentation.Arrangement

@testset "Pipeline Arrangement 1" begin
	V = [0.0 1.0; 1.0 1.0; 3.0 1.0; 4.0 1.0; 2.0 2.0; 2.0 0.0];
	EV = SparseArrays.sparse(Array{Int8, 2}([
		[1 0 0 1 0 0] #1->1,4
		[0 1 1 0 0 0] #2->2,3
		[0 0 0 0 1 1] #3->5,6
		[0 1 0 0 1 0] #4->2,5
		[0 0 1 0 0 1] #5->3,6
	]))
	bigPI = Lar.spaceindex((convert(Lar.Points, V'), Lar.cop2lar(EV)));
	
	@testset "intersect_edges" begin
		@test Lara.intersect_edges(V, EV[1, :], EV[3, :]) == [([2.0 1.0], 0.5)];
		@test Lara.intersect_edges(V, EV[1, :], EV[3, :]) == Lara.intersect_edges(V, EV[2, :], EV[3, :])
		@test Lara.intersect_edges(V, EV[1, :], EV[2, :]) == [([1.0 1.0], 0.25); ([3.0 1.0], 0.75)];
		@test Lara.intersect_edges(V, EV[2, :], EV[1, :]) == Lara.intersect_edges(V, EV[4, :], EV[5, :]);
		@test Lara.intersect_edges(V, EV[1, :], EV[4, :]) == [([1.0 1.0], 0.25)];
		@test Lara.intersect_edges(V, EV[2, :], EV[4, :]) == [([1.0 1.0], 0.0)];

	end

	@testset "merge_vertices" begin
		V1, EV1 = Lara.frag_edge(V, EV, 1, bigPI);
		ans1 = SparseArrays.sparse(Array{Int8, 2}([
				[1 0 1 0 0] #1->1,3
				[0 0 1 0 1] #2->3,5
				[0 0 0 1 1] #3->4,5
				[0 1 0 1 0] #4->2,4
			]));
		ans2 = SparseArrays.sparse(Array{Int8, 2}([
				[1 0 1] #1->1,3
				[0 1 1] #2->2,3
			]));
		@test Lara.merge_vertices!(V1, EV1)[1] == [0.0 1.0; 4.0 1.0; 1.0 1.0; 3.0 1.0; 2.0 1.0];
		@test Lara.merge_vertices!(V1, EV1)[2] == ans1;

		@test Lara.merge_vertices!(V1, EV1, [[-1]], 1)[1] == [0.0 1.0; 4.0 1.0; 2.0 1.0];
		@test Lara.merge_vertices!(V1, EV1, [[-1]], 1)[2] == ans2;
	end

	@testset "Complete Pipeline 1" begin
		Vans = [0.0 1.0; 4.0 1.0; 1.0 1.0; 3.0 1.0; 2.0 1.0; 2.0 2.0; 2.0 0.0];
		EVans = SparseArrays.sparse(Array{Int8, 2}([
				[1 0 1 0 0 0 0] # 1->3
				[0 0 1 0 1 0 0] # 3->4
				[0 0 0 1 1 0 0] # 4->5
				[0 1 0 1 0 0 0] # 2->4
				[0 0 0 0 1 1 0] # 5->6
				[0 0 0 0 1 0 1] # 5->7
				[0 0 1 0 0 1 0] # 3->6
				[0 0 0 1 0 0 1] # 4->7
			]));
		@test Lara.planar_arrangement_1(V, EV)[1] == Vans
		@test Lara.planar_arrangement_1(V, EV)[2] == EVans
	end
end

#--------------------------------------------------------------------------------------------------------------------------------

@testset "Biconnected Components" begin
	@testset "Normal Execution" begin
		EV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0 0 0 0 0 0 0 0] #1 -> 1,2		|
			[0 1 1 0 0 0 0 0 0 0] #2 -> 2,3		|
			[1 0 1 0 0 0 0 0 0 0] #3 -> 3,1		|
			[0 0 1 1 0 0 0 0 0 0] #4 -> 3,4		 |
			[0 0 1 0 0 1 0 0 0 0] #5 -> 3,6		 |
			[0 0 1 0 0 0 1 0 0 0] #6 -> 3,7		-
			[0 0 0 1 1 0 0 0 0 0] #7 -> 4,5		 |
			[0 0 0 0 1 1 0 0 0 0] #8 -> 5,6		 |
			[0 0 0 0 0 0 0 1 1 0] #9 -> 8,9		  |
			[0 0 0 0 0 0 0 0 1 1] #10-> 9,10	  |
			[0 0 0 0 0 0 0 1 0 1] #11-> 8,10	  |
		]));
		bicomp = sort(Lara.biconnected_components(EV));
		
		@test length(Lara.biconnected_components(EV)) == 3;
		@test sort(bicomp[1]) == [1, 2, 3];
		@test sort(bicomp[2]) == [4, 5, 7, 8];
		@test sort(bicomp[3]) == [9, 10, 11];
	end

	@testset "Void Components" begin
		EV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0 0] #1 -> 1,2
			[1 0 1 0] #2 -> 1,3
			[1 0 0 1] #3 -> 1,4
		]));

		@test length(Lara.biconnected_components(EV)) == 0;
	end

	@testset "Couple of Components" begin
		EV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0] #1 -> 1,2  |
			[1 1 0] #2 -> 1,2  |
			[1 0 1] #3 -> 1,2
		]));

		@test length(Lara.biconnected_components(EV)) == 1;
		@test sort(Lara.biconnected_components(EV)[1]) == [1, 2];
	end
end

#--------------------------------------------------------------------------------------------------------------------------------

@testset "Pipeline Arrangement 2" begin
	
	@testset "Containment Boxes on Nested Triangular" begin
		W = [0.0 0.0; 8.0 0.0; 4.0 6.0; 2.0 2.0; 6.0 2.0; 4.0 5.0; 9.0 0.0; 13.0 0.0; 11.0 3.0];
		copEV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0 0 0 0 0 0 0] #1 -> 1,2
			[0 1 1 0 0 0 0 0 0] #2 -> 2,3
			[1 0 1 0 0 0 0 0 0] #3 -> 3,1
			[0 0 0 1 1 0 0 0 0] #4 -> 4,5
			[0 0 0 0 1 1 0 0 0] #5 -> 5,6
			[0 0 0 1 0 1 0 0 0] #6 -> 4,6
			[0 0 0 0 0 0 1 1 0] #7 -> 7,8
			[0 0 0 0 0 0 0 1 1] #8 -> 8,9
			[0 0 0 0 0 0 1 0 1] #9 -> 7,9
		]));

		EVs = [
			SparseArrays.sparse(Array{Int8, 2}([
				[-1  1 0 0 0 0 0 0 0] #1 -> 1,2
				[ 0 -1 1 0 0 0 0 0 0] #2 -> 2,3
				[-1  0 1 0 0 0 0 0 0] #3 -> 3,1
			])),
			SparseArrays.sparse(Array{Int8, 2}([
				[0 0 0 -1  1 0 0 0 0] #4 -> 4,5
				[0 0 0  0 -1 1 0 0 0] #5 -> 5,6
				[0 0 0 -1  0 1 0 0 0] #6 -> 4,6
			])),
			SparseArrays.sparse(Array{Int8, 2}([
				[0 0 0 0 0 0 -1  1 0] #7 -> 7,8
				[0 0 0 0 0 0  0 -1 1] #8 -> 8,9
				[0 0 0 0 0 0 -1  0 1] #9 -> 7,9
			]))
		];

		shells = [
			SparseArrays.sparse(Array{Int8}([-1; -1;  1])),
			SparseArrays.sparse(Array{Int8}([-1; -1;  1])),
			SparseArrays.sparse(Array{Int8}([-1; -1;  1])),
		];

		bicon_comps = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
		shell_bboxes = [([0.0 0.0], [8.0 6.0]), ([2.0 2.0], [6.0 5.0]), ([6.0 4.0], [8.0 6.0])];
		precontainment = sparse(Array{Int8,2}([[0 0 0]; [1 0 0]; [1 0 0]]));
		containment_graph = sparse(Array{Int8,2}([[0 0 0]; [1 0 0]; [0 0 0]]));
		transitive = containment_graph;

		@test Lara.pre_containment_test(shell_bboxes) == precontainment;
		@test Lara.prune_containment_graph(3, W, EVs, shells, containment_graph) == containment_graph;

		Lar.Arrangement.transitive_reduction!(containment_graph)

		@test containment_graph == transitive;
	end

	@testset "Pipeline Arrangement 2 over K^4" begin
		
		# K^4 + edge
		W = [0.0 0.0; 4.0 0.0; 2.0 3.0; 2.0 1.5];
		copEV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0 0] #1 -> 1,2
			[0 1 1 0] #2 -> 2,3
			[1 0 1 0] #3 -> 3,1
			[1 0 0 1] #4 -> 1,2
			[0 1 0 1] #5 -> 2,3
			[0 0 1 1] #6 -> 3,1
		]));

		copFE = SparseArrays.sparse(Array{Int8, 2}([
			[1 0 0 1 1 0] #1 -> 1,4,5
			[0 1 0 0 1 1] #2 -> 2,5,6
			[0 0 1 1 0 1] #3 -> 3,4,6
			[1 1 1 0 0 0] #4 -> 1,2,3 External
		]));

		bicon_comps = [[1, 2, 3, 4, 5, 6]];

		@testset "External Cycle" begin
			@test Lara.get_external_cycle(W, copEV, copFE) == 4;
			@test typeof(Lara.get_external_cycle(W, copEV, copFE[1:3,:])) == Nothing
		end

		@testset "Biconnected Components" begin
			@test length(Lara.biconnected_components(copEV)) == 1;
			@test sort(Lara.biconnected_components(copEV)[1]) == bicon_comps[1];
		end

		@testset "Complete Component Graph" begin
			EVs = [SparseArrays.sparse(Array{Int8, 2}([
				[-1  1  0  0] #1 -> 1,2
				[ 0 -1  1  0] #2 -> 2,3
				[-1  0  1  0] #3 -> 3,1
				[-1  0  0  1] #4 -> 1,2
				[ 0 -1  0  1] #5 -> 2,3
				[ 0  0 -1  1] #6 -> 3,1
			]))];

			boundaries = [SparseArrays.sparse(Array{Int8, 2}([
				[ 1  0  0 -1  1  0] #1 -> 1,4,5
				[ 0  1  0  0 -1  1] #2 -> 2,5,6
				[ 0  0 -1  1  0 -1] #3 -> 3,4,6
			]))];

			@test Lara.componentgraph(W, copEV, bicon_comps)[1] == 1;
			@test size(Lara.componentgraph(W, copEV, bicon_comps)[2]) == (1, 1);
			@test Lara.componentgraph(W, copEV, bicon_comps)[3] == W;
			@test Lara.componentgraph(W, copEV, bicon_comps)[4] == EVs;
			@test Lara.componentgraph(W, copEV, bicon_comps)[5] == boundaries;
			@test Lara.componentgraph(W, copEV, bicon_comps)[6] == [sparse(Array{Int8}([-1; -1; 1; 0; 0; 0]))];
			@test Lara.componentgraph(W, copEV, bicon_comps)[7] == [([0.0 0.0], [4.0 3.0])];
		end
	end

	@testset "Pipeline Arrangement 2 over planar_tile" begin
		W = [
			0.0 0.0; 0.4 0.0; 0.5 0.0; 0.6 0.0; 1.0 0.0;
			0.0 0.4; 1.0 0.4;
			0.0 0.5; 1.0 0.5;
			0.0 0.6; 1.0 0.6;
			0.0 1.0; 0.4 1.0; 0.5 1.0; 0.6 1.0; 1.0 1.0;
			0.5 0.5
		];
		copEV = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 0 0 0  0 0  0 0  0 0  0 0 0 0 0  0]	#1  ->  1  2
			[1 0 0 0 0  1 0  0 0  0 0  0 0 0 0 0  0]	#2  ->  1  6
			[0 1 0 0 0  1 0  0 0  0 0  0 0 0 0 0  0]	#3  ->  2  6

			[0 0 0 1 1  0 0  0 0  0 0  0 0 0 0 0  0]	#4  ->  4  5
			[0 0 0 1 0  0 1  0 0  0 0  0 0 0 0 0  0]	#5  ->  4  7
			[0 0 0 0 1  0 1  0 0  0 0  0 0 0 0 0  0]	#6  ->  5  7

			[0 0 0 0 0  0 0  0 0  1 0  1 0 0 0 0  0]	#7  -> 10 12
			[0 0 0 0 0  0 0  0 0  1 0  0 1 0 0 0  0]	#8  -> 10 13
			[0 0 0 0 0  0 0  0 0  0 0  1 1 0 0 0  0]	#9  -> 12 13

			[0 0 0 0 0  0 0  0 0  0 1  0 0 0 1 0  0]	#10 -> 11 15
			[0 0 0 0 0  0 0  0 0  0 1  0 0 0 0 1  0]	#11 -> 11 16
			[0 0 0 0 0  0 0  0 0  0 0  0 0 0 1 1  0]	#12 -> 15 16

			[0 0 1 0 0  0 0  1 0  0 0  0 0 0 0 0  0]	#13 ->  3  8
			[0 0 1 0 0  0 0  0 1  0 0  0 0 0 0 0  0]	#14 ->  3  9
			[0 0 1 0 0  0 0  0 0  0 0  0 0 0 0 0  1]	#15 ->  3 17
			[0 0 0 0 0  0 0  1 0  0 0  0 0 0 0 0  1]	#16 ->  8 17
			[0 0 0 0 0  0 0  1 0  0 0  0 0 1 0 0  0]	#17 ->  8 14
			[0 0 0 0 0  0 0  0 1  0 0  0 0 0 0 0  1]	#18 ->  9 17
			[0 0 0 0 0  0 0  0 1  0 0  0 0 1 0 0  0]	#19 ->  9 14
			[0 0 0 0 0  0 0  0 0  0 0  0 0 1 0 0  1]	#20 -> 14 17
		]));

		copFE = SparseArrays.sparse(Array{Int8, 2}([
			[1 1 1  0 0 0  0 0 0  0 0 0  0 0 0 0 0 0 0 0]	#1 -> 1,2,3
			[0 0 0  1 1 1  0 0 0  0 0 0  0 0 0 0 0 0 0 0]	#2 -> 4,5,6
			[0 0 0  0 0 0  1 1 1  0 0 0  0 0 0 0 0 0 0 0]	#3 -> 7,8,9
			[0 0 0  0 0 0  0 0 0  1 1 1  0 0 0 0 0 0 0 0]	#4 -> 10,11,12
			[0 0 0  0 0 0  0 0 0  0 0 0  1 0 1 1 0 0 0 0]	#5 -> 13,15,16
			[0 0 0  0 0 0  0 0 0  0 0 0  0 1 1 0 0 1 0 0]	#6 -> 14,15,18
			[0 0 0  0 0 0  0 0 0  0 0 0  0 0 0 1 1 0 0 1]	#7 -> 16,17,20
			[0 0 0  0 0 0  0 0 0  0 0 0  0 0 0 0 0 1 1 1]	#8 -> 18,19,20
			[1 1 1  1 1 1  1 1 1  1 1 1  1 1 0 0 1 0 1 0]	#9 -> external
		]));

		bicon_comps = [
			[ 1,  2,  3],
			[ 4,  5,  6],
			[ 7,  8,  9],
			[10, 11, 12],
			[13, 14, 15, 16, 17, 18, 19, 20]
		];

		bboxes = [
			([0.0 0.0], [0.4 0.4])
			([0.6 0.0], [1.0 0.4])
			([0.0 0.6], [0.4 1.0])
			([0.6 0.6], [1.0 1.0])
			([0.0 0.0], [1.0 1.0])
		];

		@testset "External Cycle" begin
			@test Lara.get_external_cycle(W, copEV, copFE) == 9;
			@test typeof(Lara.get_external_cycle(W, copEV, copFE[1:8,:])) == Nothing
		end

		@testset "Biconnected Components" begin
			@test length(Lara.biconnected_components(copEV)) == 5;
			@test sort(sort(Lara.biconnected_components(copEV))[1]) == bicon_comps[1];
			@test sort(sort(Lara.biconnected_components(copEV))[2]) == bicon_comps[2];
			@test sort(sort(Lara.biconnected_components(copEV))[3]) == bicon_comps[3];
			@test sort(sort(Lara.biconnected_components(copEV))[4]) == bicon_comps[4];
			@test sort(sort(Lara.biconnected_components(copEV))[5]) == bicon_comps[5];
		end

		@testset "Containment Boxes" begin
			pre_containment = SparseArrays.sparse(Array{Int8, 2}([
				[0 0 0 0 1]
				[0 0 0 0 1]
				[0 0 0 0 1]
				[0 0 0 0 1]
				[0 0 0 0 0]
			]));
			@test Lara.pre_containment_test(bboxes) == pre_containment;
		end

		@testset "Complete Component Graph" begin
			shells = [
				sparse(Array{Int8}([-1;  1; -1])),
				sparse(Array{Int8}([-1;  1; -1])),
				sparse(Array{Int8}([ 1; -1;  1])),
				sparse(Array{Int8}([ 1; -1;  1])),
				sparse(Array{Int8}([ 1; -1;  0;  0;  1;  0; -1; 0]))
			];

			@test Lara.componentgraph(W, copEV, bicon_comps)[1] == 5;
			@test Lara.componentgraph(W, copEV, bicon_comps)[2] == SparseArrays.spzeros(Int8, 5, 5);
			@test Lara.componentgraph(W, copEV, bicon_comps)[3] == W;
			# @test abs(Lara.componentgraph(W, copEV, bicon_comps)[4]) == EV;
			# @test abs(Lara.componentgraph(W, copEV, bicon_comps)[5]) == boundaries;
			@test Lara.componentgraph(W, copEV, bicon_comps)[6] == shells;
			@test Lara.componentgraph(W, copEV, bicon_comps)[7] == bboxes;
		end
		

	end
end