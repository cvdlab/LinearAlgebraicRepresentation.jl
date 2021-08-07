using Test
using SparseArrays

@testset "Model Generation" begin
    G = [0.0 0.0 1.0 1.0;  0.0 1.0 0.0 1.0]
    T1 = SparseArrays.sparse(Array{Int8,2}([1 1 0 0;  0 1 1 0;  0 0 1 1]))
    T2 = SparseArrays.sparse(Array{Int8,2}([1 1 1]))

    @test isnothing(CAGD.Model())

    m = CAGD.Model(G)
    @test size(m, 0, 1) == 2
    @test size(m, 0, 2) == 4
    @test size(m, 1, 1) == 0
    @test size(m, 1, 2) == 4
    @test size(m, 2, 1) == 0
    @test size(m, 2, 2) == 0
    @test_throws BoundsError size(m, 3, 1)

    m = CAGD.Model(G, [T1, T2])
    @test size(m, 0, 1) == 2
    @test size(m, 0, 2) == 4
    @test size(m, 1, 1) == 3
    @test size(m, 1, 2) == 4
    @test size(m, 2, 1) == 1
    @test size(m, 2, 2) == 3
    @test_throws BoundsError size(m, 3, 1)

    @test_throws ArgumentError CAGD.Model(G, [T1])
    @test_throws ArgumentError CAGD.Model(G, [T1, T1])
    @test_throws ArgumentError CAGD.Model(G, [T1, T1, T1])

end

@testset "Model Cloning" begin
    G = [0.0 0.0 1.0 1.0;  0.0 1.0 0.0 1.0]
    T1 = SparseArrays.sparse(Array{Int8,2}([1 1 0 0;  0 1 1 0;  0 0 1 1]))
    T2 = SparseArrays.sparse(Array{Int8,2}([1 1 1]))
    m1 = CAGD.Model(G, [T1, T2])
    m2 = deepcopy(m1)
    m3 = copy(m2)
    CAGD.addModelCell!(m2, 1, [1, 3])

    @test m1.T[1] != m2.T[1]
    @test m2.T[1] == m3.T[1]
end

@testset "Topology Manipulation" begin

    G = [
        0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
        0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0
        0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0
    ]
    m = CAGD.Model(G)

    @testset "add Model Cells" begin
        CAGD.addModelCell!(m, 1, SparseArrays.sparse(Int8[-1; 1; 0; 0;  0; 0; 0; 0]))
        @test_throws MethodError CAGD.addModelCell!(m, 1, SparseArrays.sparse(Int8[-1 0 1 0  0 0 0 0]))
        @test_throws ArgumentError CAGD.addModelCell!(m, 1, SparseArrays.sparse(Int8[-1; 1; 0; 0;  0; 0; 0]))

        CAGD.addModelCells!(m, 1, SparseArrays.sparse(Int8[-1 0 1 0  0 0 0 0]))
        @test_throws ArgumentError CAGD.addModelCells!(m, 1, SparseArrays.sparse(Int8[-1 0 1 0  0 0 0]))
        @test size(m, 1) == (2, 8)

        @test_throws ArgumentError CAGD.addModelCells!(m, 0, SparseArrays.sparse(Int8[-1 0 1 0  0 0 0 0]))
        @test_throws ArgumentError CAGD.addModelCells!(m, 4, SparseArrays.sparse(Int8[-1 0 1 0  0 0 0 0]))
    end

    @testset "add Model Cells (EV)" begin
        CAGD.addModelCell!(m, 1, [2; 3], signed = true)
        @test m.T[1][3, :].nzval == [-1; 1]
        @test_throws ArgumentError CAGD.addModelCell!(m, 2, [1, 2, 4])
        CAGD.addModelCells!(m, 2, [[1, 2, 3]], signed = true)
        @test Array(m.T[2][1, :]) == [1; -1; 1]
        CAGD.addModelCells!(m, 1, [[2, 4], [3, 4]], signed = true)
        @test size(m, 2) == (1, 5)
        CAGD.addModelCell!(m, 2, [1; 2; 4; 5])
        @test m.T[2][2, :].nzval == [1; 1; 1; 1]
        CAGD.addModelCells!(m, 2, [[3; 4; 5]], signed = true)
        @test size(m, 3) == (0, 3)
    end

    @testset "get Model Lo Cell" begin
        @test CAGD.getModelLoCell(m, 2, 1) == [1; 2; 3]
        @test CAGD.getModelLoCell(m, 2, 2) == [1; 2; 4; 5]
        @test CAGD.getModelLoCell(m, 1, 2) == [1; 3]
        @test_throws ArgumentError CAGD.getModelLoCell(m, 0, 2)
        @test_throws ArgumentError CAGD.getModelLoCell(m, 4, 2)
        @test_throws ArgumentError CAGD.getModelLoCell(m, 3, 2)
        @test_throws ArgumentError CAGD.getModelLoCell(m, 2, -1)

        @test CAGD.getModelLoCell(m, 2, [1;2]) == [1; 2; 3; 4; 5]
        @test CAGD.getModelLoCell(m, 1, [1;2]) == [1; 2; 3]
        @test_throws ArgumentError CAGD.getModelLoCell(m, 0, [1; 2])
        @test_throws ArgumentError CAGD.getModelLoCell(m, 4, [1; 2])
        @test_throws ArgumentError CAGD.getModelLoCell(m, 2, [6; 2])
        @test_throws ArgumentError CAGD.getModelLoCell(m, 2, [2; 6])
        @test_throws ArgumentError CAGD.getModelLoCell(m, 2, [-1; 1])
    end

    @testset "get Model Cell Vertices" begin
        @test CAGD.getModelCellVertices(m, 2, 3) == [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0]
        ]

        @test CAGD.getModelCellVertices(m, 2, [1; 3]) == [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 1.0]
        ]

        @test CAGD.getModelCellVertices(m, 2, [3], true)[2] == [2; 3; 4]

        @test_throws ArgumentError CAGD.getModelCellVertices(m, 0, [1; 2])
        @test_throws ArgumentError CAGD.getModelCellVertices(m, 4, [1; 2])
        @test_throws ArgumentError CAGD.getModelCellVertices(m, 2, [6; 2])
        @test_throws ArgumentError CAGD.getModelCellVertices(m, 2, [2; 6])
        @test_throws ArgumentError CAGD.getModelCellVertices(m, 2, [-1; 1])
    end

    @testset "get Model Cell Geometry" begin
        @test CAGD.getModelCellGeometry(m, 2, 3) == [
            0.0 0.0 0.0
            0.0 1.0 1.0
            1.0 0.0 1.0
        ]

        @test CAGD.getModelCellGeometry(m, 2, [1; 3]) == [
            0.0 0.0 0.0 0.0
            0.0 0.0 1.0 1.0
            0.0 1.0 0.0 1.0
        ]

        @test CAGD.getModelCellGeometry(m, 2, [3], true)[2] == [2; 3; 4]

        @test_throws ArgumentError CAGD.getModelCellGeometry(m, 0, [1; 2])
        @test_throws ArgumentError CAGD.getModelCellGeometry(m, 4, [1; 2])
        @test_throws ArgumentError CAGD.getModelCellGeometry(m, 2, [6; 2])
        @test_throws ArgumentError CAGD.getModelCellGeometry(m, 2, [2; 6])
        @test_throws ArgumentError CAGD.getModelCellGeometry(m, 2, [-1; 1])
    end

    @testset "delete Model Cells" begin
        CAGD.deleteModelCell!(m, 2, 2)
        @test size(m, 2) == (2, 5)
        CAGD.addModelCell!(m, 2, [1; 2; 4; 5])
        CAGD.deleteModelCells!(m, 1, [4; 5])
        @test size(m, 2) == (1, 3)
        @test_throws ArgumentError CAGD.deleteModelCells!(m, 1, [4; 5])
        @test_throws ArgumentError CAGD.deleteModelCells!(m, 0, [4; 5])
        @test_throws ArgumentError CAGD.deleteModelCells!(m, 4, [4; 5])
    end
    
end

@testset "Geometry Manipulation" begin
    G1 = [0. 1.; 0. 1.];
    G2 = [1. 0.; 0. 1.];
    m = CAGD.Model(G1);

    @testset "add model vertices" begin
        CAGD.addModelCell!(m, 1, [1, 2]);
        CAGD.addModelVertices!(m, G2)
        @test m.G == [G1 G2]
        @test size(m, 1) == (1, 4)
        @test m.T[1] == SparseArrays.sparse(Int8[1 1 0 0])
    end

    @testset "delete model vertices" begin
        CAGD.deleteModelVertex!(m, 3)
        @test m.G == [0. 1. 0.; 0. 1. 1.];
        @test size(m, 1) == (1, 3)

        CAGD.deleteModelVertices!(m, [1, 3]);
        @test m.G == [1.; 1.][:, :]
        @test size(m, 1) == (0, 1)

        @test_throws ArgumentError CAGD.deleteModelVertices!(m, [1, 3])
        @test_throws ArgumentError CAGD.deleteModelVertices!(m, [0, 1])

    end
end

@testset "Model Manipulation" begin
    
    @testset "Purge Model" begin
        m = CAGD.Model([0.0 0.0 1.0 0.0 2.0; 0.0 1.0 0.0 2.0 0.0], [
            SparseArrays.sparse(Int8[-1 1 0 0 0; -1 0 1 0 0; 0 -1 1 0 0; 0 -1 0 1 0]),
            SparseArrays.sparse(Int8[1 -1 1 0 ])
        ])
        CAGD.modelPurge!(m)
        @test size(m, 1) == (4, 4)

        CAGD.modelPurge!(m, depth = 1)
        @test size(m, 1) == (3, 3)
    end

    @testset "Merge Model Vertices" begin
        m = CAGD.Model([
                0.0  0.0  1.0  0.0   0.0  -1.0  0.0  -1.0  -1.0
                0.0  1.0  0.0  0.0  -1.0   0.0  1.0   0.0   1.0
            ], [
                SparseArrays.sparse(Int8[
                    -1   1  0   0   0  0   0   0  0
                    -1   0  1   0   0  0   0   0  0
                     0  -1  1   0   0  0   0   0  0
                     0   0  0  -1   1  0   0   0  0
                     0   0  0  -1   0  1   0   0  0
                     0   0  0   0  -1  1   0   0  0
                     0   0  0   0   0  0  -1   1  0
                     0   0  0   0   0  0  -1   0  1
                     0   0  0   0   0  0   0  -1  1
                ]),
                SparseArrays.sparse(Int8[
                    1  -1  1  0   0  0  0   0  0
                    0   0  0  1  -1  1  0   0  0
                    0   0  0  0   0  0  1  -1  1
                ]),
        ])

        m1 = CAGD.mergeModelVertices(m, signed_merge = false);
        m2 = CAGD.mergeModelVertices(m, signed_merge = true);

        @test size(m1, 1) == (9, 6)
        @test size(m2, 1) == (9, 6)
        @test size(m1, 2) == (3, 9)
        @test size(m2, 2) == (3, 9)

        @test m1.G == m2.G
        @test sortslices(m1.T[1], dims = 1) == sortslices(abs.(m2.T[1]), dims = 1)

        @test m1.G == [
            0.0  0.0  1.0   0.0  -1.0  -1.0
            0.0  1.0  0.0  -1.0   0.0   1.0
        ]

        @test Matrix(m2.T[1]) == [
            -1   0  0   0   1  0
            -1   0  0   1   0  0
            -1   0  1   0   0  0
            -1   1  0   0   0  0
             0  -1  0   0   0  1
             0  -1  0   0   1  0
             0  -1  1   0   0  0
             0   0  0  -1   1  0
             0   0  0   0  -1  1
        ]

        @test Matrix(m2.T[2]) == [
            1  -1  0   0  0   0   0  -1   0
            0   0  1  -1  0   0  -1   0   0
            0   0  0   0  1  -1   0   0  -1
        ]
    end
end

@testset "Models Interaction" begin
    @testset "Unite Models" begin
        m1 = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
            SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
            SparseArrays.sparse(Int8[1 -1 1])
        ])
        
        m2 = CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
            SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
            SparseArrays.sparse(Int8[1 -1 1])
        ])

        m3 = CAGD.Model([0.0 -1.0 -1.0; 1.0 0.0 1.0], [
            SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
            SparseArrays.sparse(Int8[1 -1 1])
        ])
        
        m4 = CAGD.uniteModels(m1, m2)
        CAGD.uniteModels!(m1, m2)

        @test m1 == m4
        @test size(m1, 1) == (6, 6)
        @test size(m1, 2) == (2, 6)
        
        models = [
            CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ]);

            CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ]);
            
            CAGD.Model([0.0 -1.0 -1.0; 1.0 0.0 1.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ])
        ];

        CAGD.uniteModels!(m1, m3)

        @test m1 == CAGD.uniteMultipleModels(models)
    end

    @testset "Merge Models" begin
        models = [
            CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ]);
        
            CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ]);
            
            CAGD.Model([0.0 -1.0 -1.0; 1.0 0.0 1.0], [
                SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
                SparseArrays.sparse(Int8[1 -1 1])
            ])
        ]

        m = CAGD.mergeMultipleModels(models, signed_merge = true)

        @test size(m, 1) == (9, 6)
        @test size(m, 2) == (3, 9)

        @test m.G == [
            0.0  0.0  1.0   0.0  -1.0  -1.0
            0.0  1.0  0.0  -1.0   0.0   1.0
        ]

        @test Matrix(m.T[1]) == [
            -1   0  0   0   1  0
            -1   0  0   1   0  0
            -1   0  1   0   0  0
            -1   1  0   0   0  0
             0  -1  0   0   0  1
             0  -1  0   0   1  0
             0  -1  1   0   0  0
             0   0  0  -1   1  0
             0   0  0   0  -1  1
        ]

        @test Matrix(m.T[2]) == [
            1  -1  0   0  0   0   0  -1   0
            0   0  1  -1  0   0  -1   0   0
            0   0  0   0  1  -1   0   0  -1
        ]
    end
end

@testset "Model Export" begin
    modelOut = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ]);
    
    #CAGD.jlexportModel(modelOut, "./testModelExport.jl", "modelIn");

    #include("./testModelExport.jl");

    #@test modelIn == modelOut

    #rm("./testModelExport.jl")
end