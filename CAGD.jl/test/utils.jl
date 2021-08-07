using Test
CAGD = CAGD


@testset "Edge Ordering" begin
    model = CAGD.Model([
        0.0  1.0  1.0  0.0 -1.0 -1.0 -1.0  0.0  1.0
        0.0  0.0  1.0  1.0  1.0  0.0 -1.0 -1.0 -1.0
    ]);
    CAGD.addModelCells!(model, 1, [
        [1,2], [1,3], [1,4], [1,5], [1,6], [1,7], [1,8]
    ]);
    
    n = size(model, 1, 1);

    ord = CAGD.eval_ord_edges(model, 1);
    ONE = findfirst(isequal(1), ord);

    @test 1:n == vcat([ord[ONE:length(ord)]..., ord[1:ONE-1]]...)
end

@testset "Face Ordering" begin
    EV = [
        [1, 2],
        [1, 3], [1, 4], [1, 5], [1, 6], [1, 7], [1, 8], [1, 9], [1, 10],
        [2, 3], [2, 4], [2, 5], [2, 6], [2, 7], [2, 8], [2, 9], [2, 10]
    ];

    FE = [
        [1, 2, 10], [1, 3, 11], [1, 4, 12],
        [1, 5, 13], [1, 6, 14], [1, 7, 15],
        [1, 8, 16], [1, 9, 17]
    ];

    n = length(FE);

    @testset "Face Ordering x-axis" begin
        model = CAGD.Model([
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  1.0  0.0 -1.0 -1.0 -1.0  0.0  1.0
            0.0  0.0  0.0  1.0  1.0  1.0  0.0 -1.0 -1.0 -1.0
        ])

        CAGD.addModelCells!(model, 1, EV)
        CAGD.addModelCells!(model, 2, FE)

        # 1 -> n
        ord = CAGD.eval_ord_faces(model, 1)
        ONE = findfirst(isequal(1), ord)
        @test 1:n == vcat([ord[ONE:length(ord)]..., ord[1:ONE-1]]...)

        # n -> 1
        model.G[:,1:2] = model.G[:,2:-1:1]
        ord = CAGD.eval_ord_faces(model, 1)
        @test n:-1:1 == vcat([ord[ONE+1:length(ord)]..., ord[1:ONE]]...)
    end

    @testset "Face Ordering y-axis" begin
        model = CAGD.Model([
            0.0  0.0 +1.0 +1.0  0.0 -1.0 -1.0 -1.0  0.0 +1.0
            0.0 +1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  0.0 -1.0 -1.0 -1.0  0.0 +1.0 +1.0 +1.0
        ])

        CAGD.addModelCells!(model, 1, EV)
        CAGD.addModelCells!(model, 2, FE)

        # 1 -> n
        ord = CAGD.eval_ord_faces(model, 1)
        ONE = findfirst(isequal(1), ord)
        @test 1:n == vcat([ord[ONE:length(ord)]..., ord[1:ONE-1]]...)

        # n -> 1
        model.G[:,1:2] = model.G[:,2:-1:1]
        ord = CAGD.eval_ord_faces(model, 1)
        @test n:-1:1 == vcat([ord[ONE+1:length(ord)]..., ord[1:ONE]]...)
    end

    @testset "Face Ordering z-axis" begin
        model = CAGD.Model([
            0.0  0.0  1.0  1.0  0.0 -1.0 -1.0 -1.0  0.0  1.0
            0.0  0.0  0.0  1.0  1.0  1.0  0.0 -1.0 -1.0 -1.0
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
        ])

        CAGD.addModelCells!(model, 1, EV)
        CAGD.addModelCells!(model, 2, FE)

        # 1 -> n
        ord = CAGD.eval_ord_faces(model, 1)
        ONE = findfirst(isequal(1), ord)
        @test 1:n == vcat([ord[ONE:length(ord)]..., ord[1:ONE-1]]...)

        # n -> 1
        model.G[:,1:2] = model.G[:,2:-1:1]
        ord = CAGD.eval_ord_faces(model, 1)
        @test n:-1:1 == vcat([ord[ONE+1:length(ord)]..., ord[1:ONE]]...)
    end

    @testset "General Orientation" begin
        V = [
            0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
            0.0  0.0  1.0  1.0  0.0 -1.0 -1.0 -1.0  0.0  1.0
            0.0  0.0  0.0  1.0  1.0  1.0  0.0 -1.0 -1.0 -1.0
        ]
        # Generate a random basis
        T = CAGD.build_projection_matrix(rand(3,3));

        model = CAGD.Model((T*[V; ones(1, 10)])[1:3, :]);

        CAGD.addModelCells!(model, 1, EV)
        CAGD.addModelCells!(model, 2, FE)
        
        # 1 -> n
        ord = CAGD.eval_ord_faces(model, 1)
        ONE = findfirst(isequal(1), ord)
        @test 1:n == vcat([ord[ONE:length(ord)]..., ord[1:ONE-1]]...)

        # n -> 1
        model.G[:,1:2] = model.G[:,2:-1:1]
        ord = CAGD.eval_ord_faces(model, 1)
        ONE = findfirst(isequal(1), ord)
        @test n:-1:1 == vcat([ord[ONE+1:length(ord)]..., ord[1:ONE]]...)
    end

    #TODO collinear faces (manually tester)

end