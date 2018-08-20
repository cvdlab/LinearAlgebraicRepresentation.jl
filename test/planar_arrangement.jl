using Test
using LinearAlgebraicRepresentation
using SparseArrays
using LinearAlgebra
using LinearAlgebraicRepresentation.Arrangement

@testset "Edge fragmentation tests" begin
    V = [2 2; 4 2; 3 3.5; 1 3; 5 3; 1 2; 5 2]
    EV = sparse(Array{Int8, 2}([
        [1 1 0 0 0 0 0] #1->1,2
        [0 1 1 0 0 0 0] #2->2,3
        [1 0 1 0 0 0 0] #3->1,3
        [0 0 0 1 1 0 0] #4->4,5
        [0 0 0 0 0 1 1] #5->6,7
    ]))

    @testset "intersect_edges" begin
        inters1 = LinearAlgebraicRepresentation.Arrangement.intersect_edges(V, EV[5, :], EV[1, :])
        inters2 = LinearAlgebraicRepresentation.Arrangement.intersect_edges(V, EV[1, :], EV[4, :])
        inters3 = LinearAlgebraicRepresentation.Arrangement.intersect_edges(V, EV[1, :], EV[2, :])
        @test inters1 == [([2. 2.], 1/4),([4. 2.], 3/4)]
        @test inters2 == []
        @test inters3 == [([4. 2.], 1)]
    end

    @testset "frag_edge" begin
        rV, rEV = LinearAlgebraicRepresentation.Arrangement.frag_edge(V, EV, 5)
        @test rV == [1.0 2.0; 5.0 2.0; 2.0 2.0; 4.0 2.0; 4.0 2.0; 2.0 2.0]
        @test Matrix(rEV) == [1 0 0 0 0 1;
                             0 0 0 0 1 1; 
                             0 1 0 0 1 0]
    end
end

@testset "merge_vertices test set" begin
    n0 = 1e-12
    n1l = 1-1e-12
    n1u = 1+1e-12
    V = [ n0  n0; -n0  n0;  n0 -n0; -n0 -n0;
          n0 n1u; -n0 n1u;  n0 n1l; -n0 n1l;
         n1u n1u; n1l n1u; n1u n1l; n1l n1l;
         n1u  n0; n1l  n0; n1u -n0; n1l -n0]
    EV = Int8[1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
              0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
              0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0;
              0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0;
              0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0;
              0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;
              0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0;
              0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0;
              0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0;
              0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0;
              0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0;
              0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1;
              1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
              0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
              0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0;
              0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1]
    EV = sparse(EV)
    V, EV = LinearAlgebraicRepresentation.Arrangement.merge_vertices!(V, EV, [])

    @test V == [n0 n0; n0 n1u; n1u n1u; n1u n0]
    @test Matrix(EV) == [1 1 0 0;
                       0 1 1 0;
                       0 0 1 1;
                       1 0 0 1]
end

@testset "biconnected_components test set" begin
    EV = Int8[0 0 0 1 0 0 0 0 0 0 1 0; #1
              0 0 1 0 0 1 0 0 0 0 0 0; #2
              0 0 0 0 0 0 1 0 0 1 0 0; #3
              1 0 0 0 1 0 0 0 0 0 0 0; #4
              0 0 0 1 0 0 0 1 0 0 0 0; #5
              0 0 1 0 0 0 0 0 1 0 0 0; #6
              0 1 0 0 0 0 0 0 0 1 0 0; #7
              0 0 0 0 1 0 0 0 0 0 0 1; #8
              0 0 0 0 0 0 0 1 0 0 1 0; #9
              0 0 0 0 0 1 0 0 1 0 0 0; #10
              0 1 0 0 0 0 1 0 0 0 0 0; #11
              0 0 0 0 1 0 0 0 0 0 1 0; #12
              0 0 0 0 1 0 0 0 1 0 0 0] #13
    EV = sparse(EV)

    bc = LinearAlgebraicRepresentation.Arrangement.biconnected_components(EV)
    bc = Set(map(Set, bc))

    @test bc == Set([Set([1,5,9]), Set([2,6,10]), Set([3,7,11])])
end

@testset "Face creation" begin
    @testset "External cell individuation" begin
        V = [ .5 .5;  1.5   1;  1.5  2; 
             2.5  2;  2.5   1;  3.5 .5;
             3.5  3;    2 2.5;   .5  3]
    
        EV = Int8[-1  1  0  0  0  0  0  0  0;
                   0 -1  1  0  0  0  0  0  0;
                   0  0 -1  1  0  0  0  0  0;
                   0  0  0 -1  1  0  0  0  0;
                   0  0  0  0 -1  1  0  0  0;
                   0  0  0  0  0 -1  1  0  0;
                   0  0  0  0  0  0 -1  1  0;
                   0  0  0  0  0  0  0 -1  1;
                  -1  0  0  0  0  0  0  0  1;
                   0 -1  0  0  1  0  0  0  0]
        EV = sparse(EV)
        
        FE = Int8[ 0 -1 -1 -1  0  0  0  0  0  1;
                   1  1  1  1  1  1  1  1 -1  0;
                  -1  0  0  0 -1 -1 -1 -1  1 -1]
        FE = sparse(FE)
    
        @test LinearAlgebraicRepresentation.Arrangement.get_external_cycle(V, EV, FE) == 3
    end

    @testset "Containment test" begin
        V = [  0   0;    4   0;    4   2;   2   4;  0 4;
              .5  .5;  2.5  .5;  2.5 2.5;  .5 2.5;
               1   1;  1.5   1;    1   2;
               2   1;    2   2;  1.5   2;
             3.5 3.5;    3 3.5;  3.5   3]
        EV1 = Int8[ 0  0  0  0  0  0  0  0  0 -1  1  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0 -1  1  0  0  0  0  0  0;
                    0  0  0  0  0  0  0  0  0 -1  0  1  0  0  0  0  0  0]
        EV2 = Int8[ 0  0  0  0  0  0  0  0  0  0  0  0 -1  1  0  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0  0 -1  1  0  0  0;
                    0  0  0  0  0  0  0  0  0  0  0  0 -1  0  1  0  0  0]
        EV3 = Int8[ 0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0  0;
                    0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0]
        EV4 = Int8[-1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                    0 -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0;
                    0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0;
                   -1  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0]
        EV5 = Int8[ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  1  0;
                    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  1;
                    0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  1]
        EVs = map(sparse, [EV1, EV2, EV3, EV4, EV5])
        
        shell1 = Int8[-1 -1  1];
        shell2 = Int8[-1 -1  1];
        shell3 = Int8[-1 -1 -1  1];
        shell4 = Int8[-1 -1 -1 -1  1];
        shell5 = Int8[-1 -1  1];
        shells = map(sparsevec, [shell1, shell2, shell3, shell4, shell5])
        
        shell_bboxes = []
        n = 5
        for i in 1:n
            vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
            push!(shell_bboxes, LinearAlgebraicRepresentation.bbox(V[vs_indexes, :]))
        end
    
        graph = LinearAlgebraicRepresentation.Arrangement.pre_containment_test(shell_bboxes)
        @test graph == [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]
    
        graph = LinearAlgebraicRepresentation.Arrangement.prune_containment_graph(n, V, EVs, shells, graph)
        @test graph == [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
    end

    @testset "Transitive reduction" begin
        graph = [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
        LinearAlgebraicRepresentation.Arrangement.transitive_reduction!(graph)
        @test graph == [0 0 1 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
    end

    @testset "Cell merging" begin
        graph = [0 1; 0 0]
        V = [.25 .25; .75 .25; .75 .75; .25 .75;
               0   0;   1   0;   1   1;   0   1]
        EV1 = Int8[-1  1  0  0  0  0  0  0;
                    0 -1  1  0  0  0  0  0;
                    0  0 -1  1  0  0  0  0;
                   -1  0  0  1  0  0  0  0;
                   -1  0  1  0  0  0  0  0]
        EV2 = Int8[ 0  0  0  0 -1  1  0  0;
                    0  0  0  0  0 -1  1  0;
                    0  0  0  0  0  0 -1  1;
                    0  0  0  0 -1  0  0  1]
        EVs = map(sparse, [EV1, EV2])
    
        shell1 = Int8[-1 -1 -1  1  0]
        shell2 = Int8[-1 -1 -1  1]
        shells = map(sparsevec, [shell1, shell2])
    
        boundary1 = Int8[ 1  1  0  0 -1;
                          0  0  1 -1  1]
        boundary2 = Int8[ 1  1  1 -1]
        boundaries = map(sparse, [boundary1, boundary2])
    
        shell_bboxes = []
        n = 2
        for i in 1:n
            vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
            push!(shell_bboxes, LinearAlgebraicRepresentation.bbox(V[vs_indexes, :]))
        end
    
        EV, FE = LinearAlgebraicRepresentation.Arrangement.cell_merging(2, graph, V, EVs, boundaries, shells, shell_bboxes)
    
        selector = sparse(ones(Int8, 1, 3))
    
        @test selector*FE == [0  0  0  0  0  1  1  1 -1]
    end
end