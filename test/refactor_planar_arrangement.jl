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

#--------------------------------------------------------------------------------------------------------------------------------#--------------------------------------------------------------------------------------------------------------------------------

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

@testset "Pipeline Arrangement 2 over K^4" begin
	
	# K^4 + edge
	V = [0.0 0.0; 4.0 0.0; 2.0 3.0; 2.0 1.5];
    EV = SparseArrays.sparse(Array{Int8, 2}([
		[1 1 0 0] #1 -> 1,2
		[0 1 1 0] #2 -> 2,3
		[1 0 1 0] #3 -> 3,1
		[1 0 0 1] #4 -> 1,2
		[0 1 0 1] #5 -> 2,3
		[0 0 1 1] #6 -> 3,1
	]));

	EVs = [SparseArrays.sparse(Array{Int8, 2}([
		[-1  1  0  0] #1 -> 1,2
		[ 0 -1  1  0] #2 -> 2,3
		[-1  0  1  0] #3 -> 3,1
		[-1  0  0  1] #4 -> 1,2
		[ 0 -1  0  1] #5 -> 2,3
		[ 0  0 -1  1] #6 -> 3,1
	]))];

	FE = SparseArrays.sparse(Array{Int8, 2}([
		[1 0 0 1 1 0] #1 -> 1,4,5
		[0 1 0 0 1 1] #2 -> 2,5,6
		[0 0 1 1 0 1] #3 -> 3,4,6
		[1 1 1 0 0 0] #4 -> 1,2,3 External
	]));

	boundaries = [SparseArrays.sparse(Array{Int8, 2}([
		[ 1  0  0 -1  1  0] #1 -> 1,4,5
		[ 0  1  0  0 -1  1] #2 -> 2,5,6
		[ 0  0 -1  1  0 -1] #3 -> 3,4,6
	]))];

	bicon_comps = [[1, 2, 3, 4, 5, 6]];

    @testset "get_external_cycle" begin
    	@test Lara.get_external_cycle(V, EV, FE) == 4;
		@test typeof(Lara.get_external_cycle(V, EV, FV[1:3,:])) == Nothing
	end

	@testset "Biconnected Components 2" begin
		@test length(Lara.biconnected_components(EV)) == 1;
		@test sort(Lara.biconnected_components(EV)[1]) == bicon_comps;
	end



	@test Lara.componentgraph(V, copEV, bicon_comps)[1] == 1;
	@test Lara.componentgraph(V, copEV, bicon_comps)[2] == (1, 1);
	@test Lara.componentgraph(V, copEV, bicon_comps)[3] == V;
	@test Lara.componentgraph(V, copEV, bicon_comps)[4] == EVs;
	@test Lara.componentgraph(V, copEV, bicon_comps)[5] == boundaries;
	@test Lara.componentgraph(V, copEV, bicon_comps)[6] == [sparse(Array{Int8}([-1; -1; 1; 0; 0; 0]))];
	@test Lara.componentgraph(V, copEV, bicon_comps)[7] == [([0.0 0.0], [4.0 3.0])];
	

end