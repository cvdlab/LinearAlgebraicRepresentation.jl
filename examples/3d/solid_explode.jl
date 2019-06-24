using Plasm
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Debugger, Test

function solid_explode()
    #V,(VV,EV,FV,CV) = Lar.cuboid([0.5,0.5,0.5],true,[-0.5,-0.5,-0.5])
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([3,3,3],true)
    mybox = (V,CV,FV,EV)
    #mybox = Lar.Struct([ Lar.t(0,0,1.), mybox ])

    twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), mybox ])
    #twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/3,0,0), Lar.r(0,0,pi/6), mybox ])
    V,CV,FV,EV = Lar.struct2lar(twocubes)
    Plasm.view(V,CV)

	W,bases,coboundaries = Lar.chaincomplex(V,FV,EV)
	return W,bases,coboundaries
end

V,bases,coboundaries = solid_explode()
EV,FV,CV = bases
copEV,copFE,copCF = coboundaries
W = convert(Lar.Points, V')  # W by rows

objs = Lar.lar2obj(V, [copEV,copFE,copCF])

# open("./two_cubes.obj", "w") do f
# 	write(f, objs)
# end
#
#
#
# triangulated_faces = Lar.triangulate(W, [copEV, copFE])
# FVs = convert(Array{Lar.Cells}, triangulated_faces)
# Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})
#
# EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
# model = V,EVs
# Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
#
#
#
#
#
#
# function FV2EVs(copEV::Lar.ChainOp, copFE::Lar.ChainOp)
# 	EV = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
# 	FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
# 	EVs = [[EV[e] for e in fe] for fe in FE]
# 	return EVs
# end
#
#
# julia> cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1],
# [[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]],
# [[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
#
# julia> cube_2 = Lar.Struct([Lar.t(0,0,0.5), Lar.r(0,0,pi/3), cube_1])
#
# julia> V, FV, EV = Lar.struct2lar(Lar.Struct([ cube_1, cube_2 ]))
#
# julia> V, bases, coboundaries = Lar.chaincomplex(V,FV,EV)
#
# julia> (EV, FV, CV), (copEV, copFE, copCF) = bases, coboundaries
#
# julia> FV # bases[2]
# 18-element Array{Array{Int64,1},1}:
#  [1, 3, 4, 6]
#  [2, 3, 5, 6]
#  [7, 8, 9, 10]
#  [1, 2, 3, 7, 8]
#  [4, 6, 9, 10, 11, 12]
#  [5, 6, 11, 12]
#  [1, 4, 7, 9]
#  [2, 5, 11, 13]
#  [2, 8, 10, 11, 13]
#  [2, 3, 14, 15, 16]
#  [11, 12, 13, 17]
#  [11, 12, 13, 18, 19, 20]
#  [2, 3, 13, 17]
#  [2, 13, 14, 18]
#  [15, 16, 19, 20]
#  [3, 6, 12, 15, 19]
#  [3, 6, 12, 17]
#  [14, 16, 18, 20]
#
# julia> CV # bases[3]
# 3-element Array{Array{Int64,1},1}:
#  [2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20]
#  [2, 3, 5, 6, 11, 12, 13, 17]
#  [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]
#
# julia> copEV # coboundaries[1]
# 34×20 SparseMatrixCSC{Int8,Int64} with 68 stored entries: ...
#
# julia> copFE # coboundaries[2]
# 18×34 SparseMatrixCSC{Int8,Int64} with 80 stored entries: ...
#
# julia> copCF # coboundaries[3]
# 4×18 SparseMatrixCSC{Int8,Int64} with 36 stored entries: ...
#
# objs = Lar.lar2obj(V'::Lar.Points, [coboundaries...])
#
# open("./two_cubes.obj", "w") do f
# 	write(f, objs)
# end
