using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm,SparseArrays


@testset "intersect edges" begin
	V=[0 0 ; 1 1; 1/2 0; 1/2 1];
	EV = SparseArrays.sparse(Array{Int8, 2}([
		             [1 1 0 0] #1->1,2
		             [0 0 1 1] #2->3,4       
		         ]));

	W=Lar.Arrangement.intersect_edges(V, EV[1, :], EV[2, :])
	@test W[1][1]==[0.5 0.5]
	@test W[1][2]==0.5
 
end

@testset "frag_edge" begin
	V = [ 3 1; 6 8;2 -2;9 6;1 6;10 1]
    
	EV = Int8[ 1 1 0 0 0 0 ;
        	   0 0 1 1 0 0 ;
        	   0 0 0 0 1 1 ]
	EV = sparse(EV)
	edgenum = size(EV, 1)
	rV = Lar.Points(zeros(0, 2))
	rEV = SparseArrays.spzeros(Int8, 0, 0)
	model = (convert(Lar.Points,V'),Lar.cop2lar(EV))
	bigPI = Lar.spaceindex(model::Lar.LAR)
	for i in 1:edgenum
		v, ev = Lar.Arrangement.frag_edge(V, EV, i, bigPI)
		rV, rEV = Lar.skel_merge(rV, rEV, v, ev) 
	end

	@test Lar.cop2lar(rEV) == [[1,3],[2,3],[4,6],[5,6],[7,10],[9,10],[8,9]]
end


@testset "merge_vertices!" begin

	@testset "close vertices" begin
		V = [ 0.01 0.01; 0.01 -0.01; 1.01  0.99; 1.01  1.01]
		EV = Int8[1 0 1 0 ;
		          0 1 0 1 ;
		          1 0 0 1 ;
		          0 1 1 0 ]
		EV = sparse(EV)
		V, EV = Lar.Arrangement.merge_vertices!(V, EV, [],1e-1)
		@test size(V,1)==2
		@test size(EV,1)==1
	end

	@testset "double vertices" begin
		V=[0 0; 0 0; 1 1; 1 1; 2 2 ;2 2]
		EV=Int8[1 0 1 0;
		    0 1 0 1]
		EV = sparse(EV)
		V, EV = Lar.Arrangement.merge_vertices!(V, EV, [])
		@test size(V,1)==3
		@test size(EV,1)==1
	end
end


@testset "biconnected_components" begin
	@testset "Cycle graph" begin
		EV = Int8[1 1 0 0;
			0 1 1 0;
			0 0 1 1;
			1 0 0 1] 
		EV = sparse(EV)

		bc = Lar.Arrangement.biconnected_components(EV)
		@test size(bc,1)==1
	end

	@testset "Graph with 2 biconnected components" begin
		EV = Int8[1 1 0 0 0 0;
			0 1 1 0 0 0;
			1 0 1 0 0 0;
			1 0 0 0 1 0;
			0 0 0 1 1 0;
			0 0 0 1 0 1;
			0 0 0 0 1 1] ;
		EV = sparse(EV);

		bc = Lar.Arrangement.biconnected_components(EV)
		@test size(bc,1)==2
	end
end


@testset "get_external_cycle" begin
	@testset "Example 1" begin

		V = [ 0 0; 1 0;2 0;2 1;1 1;0 1]
	    
		EV = Int8[-1  1  0  0  0  0 ;
			   0 -1  1  0  0  0 ;
			   0  0 -1  1  0  0 ;
			   0  0  0 -1  1  0 ;
			   0  0  0  0 -1  1 ;
			   1  0  0  0  0 -1 ;
			   0  -1  0  0  1  0 
			   ]
		EV = sparse(EV)
		
		FE = Int8[ 1 1 1 1 1 1 0;
			   1 0 0 0 1 1 -1;
			   0 1 1 1 0 0 -1]
		FE = sparse(FE)
		@test Lar.Arrangement.get_external_cycle(V, EV, FE) == 1
	end

	@testset "Example 2" begin

		V = [ 0 0; 2 0; 2 2; 0 2; 
		      1 0; 3 0; 3 3; 1 2]
	    
		EV = Int8[-1  1  0  0  0  0  0  0;
			   0 -1  1  0  0  0  0  0;
			   0  0 -1  1  0  0  0  0;
			  -1  0  0  1  0  0  0  0;
			   0  0  0  0 -1  1  0  0;
			   0  0  0  0  0 -1  1  0;
			   0  0  0  0  0  0 -1  1;
			   0  0  0  0 -1  0  0  1;
			   ]
		EV = sparse(EV)
		
		FE = Int8[-1 -1 -1  1  0  0  0  0;
			   0  0  0  0 -1 -1 -1  1;
			   0  0  0  1 -1 -1 -1  0]
		FE = sparse(FE)
		@test Lar.Arrangement.get_external_cycle(V, EV, FE) == 3
	end

end



@testset "Containment graph" begin
	@testset "pre_containment and prune_containment" begin

		V = [ 0 0; 2 0; 0 2; 3/2 3/2; 1/2 3/2;1/2 1/2 ; 3/2 1/2]
		    
		EV = Int8[ -1  1  0  0  0  0  0;
			   0  -1  1  0  0  0  0;
			   -1  0  1  0  0  0  0]

		EV1=Int8[  0  0  0  -1  1  0  0;
			   0  0  0  0  -1  1  0;
			   0  0  0  0  0  -1  1;
			   0  0  0  -1  0  0  1
			   ]
		EVs = map(sparse, [EV, EV1])

		shell1=Int8[-1 -1 1]
		shell2=Int8[-1 -1 -1 1]
		shells = map(sparsevec, [shell1, shell2])
		     
		shell_bboxes = []
		n = 2
		for i in 1:n
		    vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
		    push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
		end
		   
		graph = Lar.Arrangement.pre_containment_test(shell_bboxes)
		@test graph == [0 0;1 0];
		   
		graph = Lar.Arrangement.prune_containment_graph(n, V, EVs, shells, graph)
		@test graph == [0 0;0 0]
	end

	@testset "transitive_reduction" begin

		graph = [0 1 1 1 ; 0 0 1 1 ; 0 0 0 1 ; 0 0 0 0 ]
	 	Lar.Arrangement.transitive_reduction!(graph)
		@test graph == [0 1 0 0 ; 0 0 1 0; 0 0 0 1 ; 0 0 0 0 ]
	end

end

@testset "component graph" begin
	
	V = [ 0 0; 2 0; 0 1; 3/2 1/2; 3/2 1/2; 1/2 3/2; 3/2 3/2]
	    
	EV = Int8[ 1  1  0  0  0  0  0;
		   0  1  1  0  0  0  0;
		   1  0  1  0  0  0  0;
		   0  0  0  1  1  0  0;
		   0  0  0  0  1  0  1;
		   0  0  0  0  0  1  1;
		   0  0  0  1  0  1  0
		   ]
	EV = sparse(EV)
		
	bc = Lar.Arrangement.biconnected_components(EV)

	n, containment_graph, V, EVs, boundaries, shells, shell_bboxes=Lar.Arrangement.componentgraph(V, EV, bc)

	@test n==2
	@test boundaries==[[1 1 -1 -1],[1 1 -1]]
	@test shell_bboxes==[([0.5 0.5], [1.5 1.5]), ([0.0 0.0], [2.0 1.0])]
end 

@testset "Planar Arrangement" begin

	V = [ 0 0; 1 0; 1 1; 0 1; 2 1]
    
	EV=[[1,2],[2,3],[3,4],[1,4],[3,5]]
	
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)

	V, copEV, copFE = Lar.Arrangement.planar_arrangement(V::Lar.Points, cop_EW::Lar.ChainOp)

	@test Lar.cop2lar(copEV)==[[1,2],[2,3],[3,4],[1,4]]
end
