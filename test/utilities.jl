using Test
using SparseArrays
using LARLIB

@testset "Bounding boxes building test" begin
    bboxV = [.56 .28; .84 .57; .35  1.0; .22  .43]
    @test LARLIB.bbox(bboxV) == ([.22 .28], [.84 1.0])
end

@testset "Bounding boxes containment test" begin
    bboxA = ([0. 0.], [1. 1.])
    bboxB = ([.5 .5], [1.5 1.5])
    bboxC = ([1. 1.], [1.25 1.25])
    bboxD = ([0 .75], [.25 1])
    bboxE = ([0 1.25], [.25 1.5])

    @test LARLIB.bbox_contains(bboxA, bboxD)
    @test LARLIB.bbox_contains(bboxB, bboxC)
    @test !LARLIB.bbox_contains(bboxA, bboxB)
    @test !LARLIB.bbox_contains(bboxA, bboxE)
end

@testset "Face area calculation test" begin
    facetestV = Float64[2 1; 1 2; 0 0; 1 1; 2 0; 0 2]
    facetestEV = spzeros(Int8, 6, 6)
    facetestEV[1, [1, 4]] = [-1, 1]; facetestEV[2, [2, 4]] = [-1, 1]
    facetestEV[3, [2, 6]] = [-1, 1]; facetestEV[4, [3, 6]] = [-1, 1]
    facetestEV[5, [3, 5]] = [-1, 1]; facetestEV[6, [1, 5]] = [-1, 1]
    facetestFE = spzeros(Int8, 2, 6)
    facetestFE[1, :] = [ 1 -1  1 -1  1 -1]
    facetestFE[2, :] = [-1  1 -1  1 -1  1]

    @test LARLIB.face_area(facetestV, facetestEV, facetestFE[1,:]) == -LARLIB.face_area(facetestV, facetestEV, facetestFE[2,:])
end

@testset "Binaryrange" begin
    @test LARLIB.binaryRange(1) == [ "0", "1" ]
    @test LARLIB.binaryRange(2) == [ "00", "01", "10", "11"]
    @test LARLIB.binaryRange(3) == [ "000", "001", "010", "011", "100", "101", "110", "111"]
end

