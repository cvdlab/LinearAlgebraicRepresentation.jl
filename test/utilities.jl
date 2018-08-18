using Base.Test
using LinearAlgebraicRepresentation

@testset "Bounding boxes building test" begin
    V = [.56 .28; .84 .57; .35  1.0; .22  .43]
    @test LinearAlgebraicRepresentation.bbox(V) == ([.22 .28], [.84 1.0])
end

@testset "Bounding boxes containment test" begin
    bboxA = ([0. 0.], [1. 1.])
    bboxB = ([.5 .5], [1.5 1.5])
    bboxC = ([1. 1.], [1.25 1.25])
    bboxD = ([0 .75], [.25 1])
    bboxE = ([0 1.25], [.25 1.5])

    @test LinearAlgebraicRepresentation.bbox_contains(bboxA, bboxD)
    @test LinearAlgebraicRepresentation.bbox_contains(bboxB, bboxC)
    @test !LinearAlgebraicRepresentation.bbox_contains(bboxA, bboxB)
    @test !LinearAlgebraicRepresentation.bbox_contains(bboxA, bboxE)
end
@testset "Face area calculation test" begin
    V = Float64[2 1; 1 2; 0 0; 1 1; 2 0; 0 2]
    EV = spzeros(Int8, 6, 6)
    EV[1, [1, 4]] = [-1, 1]; EV[2, [2, 4]] = [-1, 1]
    EV[3, [2, 6]] = [-1, 1]; EV[4, [3, 6]] = [-1, 1]
    EV[5, [3, 5]] = [-1, 1]; EV[6, [1, 5]] = [-1, 1]
    FE = spzeros(Int8, 2, 6)
    FE[1, :] = [ 1 -1  1 -1  1 -1]
    FE[2, :] = [-1  1 -1  1 -1  1]

    @test LinearAlgebraicRepresentation.face_area(V, EV, FE[1,:]) == -LinearAlgebraicRepresentation.face_area(V, EV, FE[2,:])
end

