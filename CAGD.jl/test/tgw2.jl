using SparseArrays, Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 1 METHODS
#    frag_edge_channel
#    frag_edge
#    intersect_edges
#    merge_vertices!
#-------------------------------------------------------------------------------

@testset "Edge Fragmentation" begin
    model = CAGD.Model(hcat([
        [2., 2.], [4., 2.], [3., 3.5], [1., 3.], [5., 3.], [1., 2.], [5., 2.]
    ]...))
    CAGD.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0 0 0 0 0] #1->1,2
        [0 1 1 0 0 0 0] #2->2,3
        [1 0 1 0 0 0 0] #3->1,3
        [0 0 0 1 1 0 0] #4->4,5
        [0 0 0 0 0 1 1] #5->6,7
    ])))
    bigPI = CAGD.spaceIndex(model, 1)

    @testset "intersect_edges" begin
        V = convert(Lar.Points, model.G')
        EV = model.T[1]

        inters1 = CAGD.intersect_edges(V, EV[5, :], EV[1, :])
        inters2 = CAGD.intersect_edges(V, EV[1, :], EV[4, :])
        inters3 = CAGD.intersect_edges(V, EV[1, :], EV[2, :])

        @test inters1 == [([2. 2.], 1/4),([4. 2.], 3/4)]
        @test inters2 == []
        @test inters3 == [([4. 2.], 1)]
    end

    @testset "frag_edges" begin
        V1, EV1 = CAGD.frag_edge(model, 1, bigPI)
        V4, EV4 = CAGD.frag_edge(model, 4, bigPI)
        V5, EV5 = CAGD.frag_edge(model, 5, bigPI)

        @test V1 == [2.0 2.0; 4.0 2.0; 2.0 2.0; 4.0 2.0]
        @test EV1 == SparseArrays.sparse(Array{Int8, 2}([1 1 0 0]))

        @test isapprox(V4, [1.0 3.0; 5.0 3.0; 2.666666666 3.0; 3.33333333 3.0])
        @test EV4 == SparseArrays.sparse(Array{Int8, 2}([
            [1 0 1 0]
            [0 0 1 1]
            [0 1 0 1]
        ]))

        @test V5 == [1.0 2.0; 5.0 2.0; 2.0 2.0; 2.0 2.0; 4.0 2.0; 4.0 2.0]
        @test EV5 == SparseArrays.sparse(Array{Int8, 2}([
            [1 0 0 1 0 0]
            [0 0 0 1 0 1]
            [0 1 0 0 0 1]
        ]))
    end

    @testset "frag_edge_channel" begin
                                                                                # TODO
    end
end

@testset "merge_vertices" begin
    n0 = 1e-12
    n1l = 1-1e-12
    n1u = 1+1e-12
    model = CAGD.Model([
        n0 -n0  n0 -n0  n0 -n0  n0 -n0 n1u n1l n1u n1l n1u n1l n1u n1l
        n0  n0 -n0 -n0 n1u n1u n1l n1l n1u n1u n1l n1l  n0  n0 -n0 -n0
    ])
    CAGD.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
        [1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1]
        [1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1]
    ])))
    CAGD.merge_vertices!(model)

    @test model.G == [n0 n0; n0 n1u; n1u n1u; n1u n0]'
    @test Array(model.T[1]) == [1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1]
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT 1
#-------------------------------------------------------------------------------

@testset "Planar Arrangement 1" begin
    model = CAGD.Model([
        2.0 2.0 4.0 4.0 1.5 3.0 4.5 3.0
        2.0 4.0 4.0 2.0 3.0 4.5 3.0 1.5
    ])
    CAGD.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
        1 1 0 0 0 0 0 0
        0 1 1 0 0 0 0 0
        0 0 1 1 0 0 0 0
        1 0 0 1 0 0 0 0
        0 0 0 0 1 1 0 0
        0 0 0 0 0 1 1 0
        0 0 0 0 0 0 1 1
        0 0 0 0 1 0 0 1
    ])))
    sigma = SparseArrays.sparse(Array{Int8, 1}([
        0; 0; 0; 0; 1; 1; 1; 1
    ]))
    model, sigma, edges = CAGD.planar_arrangement_1(model, sigma, true)
    @test model.G == [
        2.0 2.0 2.0 2.0 4.0 2.5 3.5 4.0 4.0 4.0 2.5 3.5 1.5 3.0 4.5 3.0
        2.0 4.0 3.5 2.5 4.0 4.0 4.0 2.0 3.5 2.5 2.0 2.0 3.0 4.5 3.0 1.5
    ]
    @test model.T[1] == abs.(Lar.coboundary_0([
        [1,4],[4,3],[3,2],[2,6],[6,7],[7,5],
        [5,9],[9,10],[10,8],[1,11],[11,12],[12,8],
        [13,3],[3,6],[6,14],[14,7],[7,9],[9,15],
        [15,10],[10,12],[12,16],[13,4],[4,11],[11,16],
    ]))
    @test sigma == 13:24
    @test edges == [
        [ 1,  2,  3], [ 4,  5,  6], [ 7,  8,  9], [10, 11, 12],
        [13, 14, 15], [16, 17, 18], [19, 20, 21], [22, 23, 24]
    ]
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - CLEAN DECOMPOSITION
#-------------------------------------------------------------------------------

@testset "Clean Decomposition" begin
    model = CAGD.Model([
        2.0 2.0 2.0 2.0 4.0 2.5 3.5 4.0 4.0 4.0 2.5 3.5 1.5 3.0 4.5 3.0
        2.0 4.0 3.5 2.5 4.0 4.0 4.0 2.0 3.5 2.5 2.0 2.0 3.0 4.5 3.0 1.5
    ])
    CAGD.addModelCells!(model, 1, abs.(Lar.coboundary_0([
        [ 1,  4], [ 4,  3], [ 3,  2], [ 2,  6], [ 6,  7], [ 7,  5],
        [ 5,  9], [ 9, 10], [10,  8], [ 1, 11], [11, 12], [12,  8],
        [13,  3], [ 3,  6], [ 6, 14], [14,  7], [ 7,  9], [ 9, 15],
        [15, 10], [10, 12], [12, 16], [13,  4], [ 4, 11], [11, 16]
    ])))
    sigma = collect(13:24)
    edges = [
        [ 1,  2,  3], [ 4,  5,  6], [ 7,  8,  9], [10, 11, 12],
        [13, 14, 15], [16, 17, 18], [19, 20, 21], [22, 23, 24]
    ]
    model, edges = CAGD.cleandecomposition(model, sigma, edges)
    @test model.G == [
        2.0 2.0 2.5 3.5 4.0 4.0 2.5 3.5 1.5 3.0 4.5 3.0
        3.5 2.5 4.0 4.0 3.5 2.5 2.0 2.0 3.0 4.5 3.0 1.5
    ]
    @test model.T[1] == abs.(Lar.coboundary_0([
        [ 2,  1], [ 3,  4],
        [ 5,  6], [ 8,  7],
        [ 9,  1], [ 1,  3], [3, 10], [10,  4], [ 4,  5], [ 5, 11],
        [11,  6], [ 6,  8], [8, 12], [ 9,  2], [ 2,  7], [ 7, 12]
    ]))
    @test edges == [
        [         1], [         2], [         3], [         4],
        [ 5,  6,  7], [ 8,  9, 10], [11, 12, 13], [14, 15, 16]
    ]
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - BICONNECTED COMPONENTS
#-------------------------------------------------------------------------------

@testset "Biconnected Components" begin
    EV = SparseArrays.sparse(Int8[
        0 0 0 1 0 0 0 0 0 0 1 0 #1
        0 0 1 0 0 1 0 0 0 0 0 0 #2
        0 0 0 0 0 0 1 0 0 1 0 0 #3
        1 0 0 0 1 0 0 0 0 0 0 0 #4
        0 0 0 1 0 0 0 1 0 0 0 0 #5
        0 0 1 0 0 0 0 0 1 0 0 0 #6
        0 1 0 0 0 0 0 0 0 1 0 0 #7
        0 0 0 0 1 0 0 0 0 0 0 1 #8
        0 0 0 0 0 0 0 1 0 0 1 0 #9
        0 0 0 0 0 1 0 0 1 0 0 0 #10
        0 1 0 0 0 0 1 0 0 0 0 0 #11
        0 0 0 0 1 0 0 0 0 0 1 0 #12
        0 0 0 0 1 0 0 0 1 0 0 0 #13
    ])

    bc = CAGD.biconnected_components(EV)
    bc = Set(map(Set, bc))

    @test bc == Set([Set([1,5,9]), Set([2,6,10]), Set([3,7,11])])
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 - COMPONENT GRAPH (TGW) METHODS
#    get_external_cycle
#    pre_containment_test
#    prune_containment_graph
#    transitive_reduction!
#-------------------------------------------------------------------------------

@testset "TGW Methods" begin
    @testset "External cell individuation" begin
        V = [ .5 .5;  1.5   1;  1.5  2;
             2.5  2;  2.5   1;  3.5 .5;
             3.5  3;    2 2.5;   .5  3]

        EV = SparseArrays.sparse(Int8[
            -1  1  0  0  0  0  0  0  0
             0 -1  1  0  0  0  0  0  0
             0  0 -1  1  0  0  0  0  0
             0  0  0 -1  1  0  0  0  0
             0  0  0  0 -1  1  0  0  0
             0  0  0  0  0 -1  1  0  0
             0  0  0  0  0  0 -1  1  0
             0  0  0  0  0  0  0 -1  1
            -1  0  0  0  0  0  0  0  1
             0 -1  0  0  1  0  0  0  0
        ])

        FE = SparseArrays.sparse(Int8[
             0 -1 -1 -1  0  0  0  0  0  1
             1  1  1  1  1  1  1  1 -1  0
            -1  0  0  0 -1 -1 -1 -1  1 -1
        ])

        @test CAGD.get_external_cycle(V, EV, FE) == 3

        model = CAGD.Model(hcat([[0.0,0.0], [4.0,0.0], [2.0,3.0], [2.0,1.5]]...))
        CAGD.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
    		[1 1 0 0] #1 -> 1,2
    		[0 1 1 0] #2 -> 2,3
    		[1 0 1 0] #3 -> 3,1
    		[1 0 0 1] #4 -> 1,2
    		[0 1 0 1] #5 -> 2,3
    		[0 0 1 1] #6 -> 3,1
    	])))
        CAGD.addModelCells!(model, 2, SparseArrays.sparse(Array{Int8, 2}([
    		[1 0 0 1 1 0] #1 -> 1,4,5
    		[0 1 0 0 1 1] #2 -> 2,5,6
    		[0 0 1 1 0 1] #3 -> 3,4,6
    		[1 1 1 0 0 0] #4 -> 1,2,3 External
    	])))

        @test CAGD.get_external_cycle(model) == 4
        CAGD.deleteModelCell!(model, 2, 4)
        @test isnothing(CAGD.get_external_cycle(model))
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

        shell1 = Int8[-1 -1  1]
        shell2 = Int8[-1 -1  1]
        shell3 = Int8[-1 -1 -1  1]
        shell4 = Int8[-1 -1 -1 -1  1]
        shell5 = Int8[-1 -1  1]
        shells = map(sparsevec, [shell1, shell2, shell3, shell4, shell5])

        shell_bboxes = []
        n = 5
        for i in 1:n
            vs_indexes = (abs.(EVs[i]')*abs.(shells[i])).nzind
            push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
        end

        graph = CAGD.pre_containment_test(shell_bboxes)
        @test graph == [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 1 0]

        graph = CAGD.prune_containment_graph(V, EVs, shells, graph)
        @test graph == [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
    end

    @testset "Transitive reduction" begin
        graph = [0 0 1 1 0; 0 0 1 1 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
        CAGD.transitive_reduction!(graph)
        @test graph == [0 0 1 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 0; 0 0 0 0 0]
    end
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2 METHODS
#    componentgraph (aka TGW)
#    cell_merging
#-------------------------------------------------------------------------------

@testset "Component Graph" begin
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
            push!(shell_bboxes, Lar.bbox(V[vs_indexes, :]))
        end

        EV, FE = CAGD.cell_merging(2, graph, V, EVs, boundaries, shells, shell_bboxes)

        selector = sparse(ones(Int8, 1, 3))

        @test selector*FE == [0  0  0  0  0  1  1  1 -1]
    end
end

@testset "Cell Merging" begin
end

#-------------------------------------------------------------------------------
#   PLANAR ARRANGEMENT PIPELINE - PART 2
#-------------------------------------------------------------------------------

@testset "Planar 2" begin
    model = CAGD.Model([
        0.0 0 2 0 0 1 4 4 3 2 4 4
        0.0 2 0 3 4 4 2 0 0 4 4 2
    ]);
    CAGD.addModelCells!(model, 1, abs.(Lar.coboundary_0([
        [ 1,  2], [ 2,  3], [ 3,  1],
        [ 4,  5], [ 5,  6], [ 6,  7], [ 7,  8], [ 8,  9], [ 9, 4],
        [10, 11], [11, 12], [12, 10]
    ])))

    bicon_comps = [
        [ 1,  2,  3],
        [ 4,  5,  6,  7,  8,  9],
        [10, 11, 12]
    ]
    model = CAGD.planar_arrangement_2(model, bicon_comps)
end

#-------------------------------------------------------------------------------
#   COMPLETE PLANAR ARRANGEMENT
#-------------------------------------------------------------------------------

@testset "Complete Pipeline" begin

    #=== TRIFORCE MODEL DEFINITION ===#
    triforce = CAGD.Model([
        0.0  4.0  8.0  4.0  2.0  6.0  4.0 -4.0 12.0  4.0  3.0  5.0  4.0  4.0
        0.0  6.0  0.0  0.0  3.0  3.0 -6.0  6.0  6.0  1.0  2.5  2.5  1.0  2.0
    ])
    CAGD.addModelCells!(triforce, 1, SparseArrays.sparse(Array{Int8, 2}([
        [1 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 1 1 0 0 0 0 0 0 0 0 0 0 0]
        [1 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 1 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 1 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 1 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 1 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 1]
    ])))
    sigma = convert(Lar.Chain, SparseArrays.sparse([
        1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0
    ]))

    #=== PLANAR ARRANGEMENT 1 ===#
    arranged_triforce = CAGD.Model([
        0.0  4.0  2.0  8.0  6.0  4.0  4.0 -4.0 12.0  4.0  3.0  5.0  4.0
        0.0  6.0  3.0  0.0  3.0  0.0 -6.0  6.0  6.0  1.0  2.5  2.5  2.0
    ])
    CAGD.addModelCells!(arranged_triforce, 1, SparseArrays.sparse(Int8[
        1 0 1 0 0 0 0 0 0 0 0 0 0
        0 1 1 0 0 0 0 0 0 0 0 0 0
        0 1 0 0 1 0 0 0 0 0 0 0 0
        0 0 0 1 1 0 0 0 0 0 0 0 0
        1 0 0 0 0 1 0 0 0 0 0 0 0
        0 0 0 1 0 1 0 0 0 0 0 0 0
        0 0 1 0 0 1 0 0 0 0 0 0 0
        0 0 1 0 1 0 0 0 0 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0
        1 0 0 0 0 0 1 0 0 0 0 0 0
        1 0 0 0 0 0 0 1 0 0 0 0 0
        0 1 0 0 0 0 0 1 0 0 0 0 0
        0 1 0 0 0 0 0 0 1 0 0 0 0
        0 0 0 1 0 0 1 0 0 0 0 0 0
        0 0 0 1 0 0 0 0 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 1 0 0
        0 0 0 0 0 0 0 0 0 0 1 1 0
        0 0 0 0 0 0 0 0 0 1 0 1 0
        0 0 0 0 0 0 0 0 0 1 0 0 1
    ]))
    arranged_sigma = collect(1:6)
    arranged_edge_map = Array{Int64,1}[
        [ 1,  2], [ 3,  4], [ 5,  6],
        [     7], [     8], [     9],
        [10, 11], [12, 13], [14, 15],
        [    16], [    17], [    18], [    19]
    ]

    a, b, c = CAGD.planar_arrangement_1(triforce, sigma, true)
    @test arranged_triforce == a
    @test arranged_sigma    == b
    @test arranged_edge_map == c

    #=== CLEAN DECOMPOSITION ===#
    cleaned_triforce = CAGD.Model([
        0.0 4.0 2.0 8.0 6.0 4.0 4.0 3.0 5.0 4.0
        0.0 6.0 3.0 0.0 3.0 0.0 1.0 2.5 2.5 2.0
    ])
    CAGD.addModelCells!(cleaned_triforce, 1, SparseArrays.sparse(Int8[
        1 0 1 0 0 0 0 0 0 0
        0 1 1 0 0 0 0 0 0 0
        0 1 0 0 1 0 0 0 0 0
        0 0 0 1 1 0 0 0 0 0
        1 0 0 0 0 1 0 0 0 0
        0 0 0 1 0 1 0 0 0 0
        0 0 1 0 0 1 0 0 0 0
        0 0 1 0 1 0 0 0 0 0
        0 0 0 0 1 1 0 0 0 0
        0 0 0 0 0 0 1 1 0 0
        0 0 0 0 0 0 0 1 1 0
        0 0 0 0 0 0 1 0 1 0
        0 0 0 0 0 0 1 0 0 1
    ]))
    cleaned_edge_map = Array{Int64,1}[
        [1, 2], [3, 4], [5, 6],
        [   7], [   8], [   9],
        [    ], [    ], [    ],
        [  10], [  11], [  12], [  13]
    ]

    a, b = CAGD.cleandecomposition(a, b, c)
    @test cleaned_triforce == a
    @test cleaned_edge_map == b

    #=== BICONNECTED COMPONENTS ===#
    bicon_comps = Array{Int64,1}[[8, 9, 7, 5, 6, 4, 3, 2, 1], [12, 11, 10]]
    c = CAGD.biconnected_components(a.T[1])
    @test bicon_comps == c

    #=== REMOVE DANGLING EDGES FROM EDGE MAP ===#
    b = CAGD.remove_mapping_dangling_edges(b, c)
    purged_edge_map = [
        [1, 2], [3, 4], [5, 6],
        [   7], [   8], [   9],
        [    ], [    ], [    ],
        [  10], [  11], [  12], [    ]
    ]
    @test purged_edge_map == b

    #=== PLANAR ARRANGEMENT 2 ===#

    rearranged_triforce = CAGD.Model([
        0.0 4.0 2.0 8.0 6.0 4.0 4.0 3.0 5.0
        0.0 6.0 3.0 0.0 3.0 0.0 1.0 2.5 2.5
    ])
    CAGD.addModelCells!(rearranged_triforce, 1, SparseArrays.sparse(Int8[
        -1  0  1  0  0  0  0  0  0
         0 -1  1  0  0  0  0  0  0
         0 -1  0  0  1  0  0  0  0
         0  0  0 -1  1  0  0  0  0
        -1  0  0  0  0  1  0  0  0
         0  0  0 -1  0  1  0  0  0
         0  0 -1  0  0  1  0  0  0
         0  0 -1  0  1  0  0  0  0
         0  0  0  0 -1  1  0  0  0
         0  0  0  0  0  0 -1  1  0
         0  0  0  0  0  0  0 -1  1
         0  0  0  0  0  0 -1  0  1
    ]))
    CAGD.addModelCells!(rearranged_triforce, 2, SparseArrays.sparse(Int8[
        0  0  0  0  0  0  1 -1 -1  1  1 -1
        -1 0  0  0  1  0 -1  0  0  0  0  0
        0  1 -1  0  0  0  0  1  0  0  0  0
        0  0  0  1  0 -1  0  0  1  0  0  0
        0  0  0  0  0  0  0  0  0 -1 -1  1
    ]))
    a = CAGD.planar_arrangement_2(a, c)
    @test rearranged_triforce == a

end

#-------------------------------------------------------------------------------
#   POSSIBLY ERRORING DATA
#-------------------------------------------------------------------------------

@testset "Formelly Erroring Data" begin
    @testset "Non-Manifold Triangles" begin
        model = CAGD.Model([
            0.0 1.0 0.0 0.5 0.1
            0.0 0.0 1.0 0.1 0.5
        ])
        CAGD.addModelCells!(model, 1, SparseArrays.sparse(Array{Int8, 2}([
            -1  1  0  0  0
            -1  0  1  0  0
             0 -1  1  0  0
            -1  0  0  1  0
            -1  0  0  0  1
             0  0  0 -1  1
        ])))
        arrangedModel = deepcopy(model)
        CAGD.addModelCells!(model, 2, SparseArrays.sparse(Array{Int8, 2}([
             1 -1  1  0  0  0
             0  0  0  1 -1  1
        ])))
        CAGD.addModelCells!(arrangedModel, 2, SparseArrays.sparse(Array{Int8,2}([
             1 -1  1 -1  1 -1
             0  0  0  1 -1  1
        ])))
        @test arrangedModel == CAGD.planar_arrangement(model)

    end
end


#-------------------------------------------------------------------------------
#   BIG RANDOM DATASETS
#-------------------------------------------------------------------------------

@testset "Planar Arrangement - Big Random Data" begin
    npts = 100
    ndgs = 200
    V = rand(2, npts)
    model = CAGD.Model(rand(2, npts))
    cv = unique([rand(1 : npts, 2) for i = 1 : 10*npts])[1 : ndgs*2]
    cv = [[c[1] c[2]] for c in cv if c[1] != c[2]][1 : ndgs]
    I = []
    J = []
    K = ones(Int8, 2*ndgs)
    for i = 1 : ndgs
        push!(J, cv[i]...)
        push!(I, [i; i]...)
    end
    CAGD.addModelCells!(model, 1, SparseArrays.sparse(I, J, K, ndgs, npts))
    model1, edge_map = CAGD.planar_arrangement_1(model, spzeros(Int8, 0), true)
    bicon_comps = CAGD.biconnected_components(model1.T[1])
    edge_map = CAGD.remove_mapping_dangling_edges(edge_map, bicon_comps)
    model2 = CAGD.planar_arrangement_2(model1, bicon_comps)

    @testset "Closure of 2-cells" begin
        for i = 1 : size(model2, 2, 1)
            edges  = findnz(model2.T[2][i, :])[1]
            points = abs.(model2.T[1][edges, :])
            @test sum(points) == 2 * length(edges)
        end
    end

    @testset "Well formed cells" begin
        @test sum(sum(abs.(model2.T[1]), dims = 1) .<= 1) == 0
        @test sum(sum(abs.(model2.T[1]), dims = 2) .!= 2) == 0
        @test sum(sum(abs.(model2.T[2]), dims = 1) .>= 3) == 0
        @test sum(sum(abs.(model2.T[2]), dims = 2) .<= 2) == 0
    end

    @testset "Euler" begin
        @test size(model2, 0, 2) - size(model2, 1, 1) + size(model2, 2, 1) == 1
    end

    @testset "1-Cell Compactness" begin
        #@test sum(abs.(model2.T[1])) == 2 * size(model, 1, 2)
    end

    @testset "2-Cell Compactness" begin
        #@test sum(abs.(model2.T[2])) == 2 * size(model2, 2, 2)
    end

    @testset "Biconnected Components" begin
        # for comp in bicon_comps
        #     sum(
        #         (sum(abs.(model.T[1][comp, :]), dims=1) .!= 0)
        #         + (sum(abs.(model.T[1][setdiff(1:size(model, 1, 1), comp), :]), dims=1) .!= 0)
        #         .!= 1
        #     ) == 0
        # end
    end
end
