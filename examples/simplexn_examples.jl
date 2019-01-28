using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

VOID = [[]], [[1]]


# example 1

V = [[0,0] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];
FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];
pattern = repeat([1,2,-3],outer=4);
model = (V,FV);
W,FW = Lar.extrudeSimplicial(model, pattern);
Plasm.view(W,FW)

# example 2

model = Lar.extrudeSimplicial( VOID, ones(10) )
Plasm.view(model)
model = Lar.extrudeSimplicial( model, ones(10) )
Plasm.view(model)
model = Lar.extrudeSimplicial( model, ones(10) )
Plasm.view(model)


# example 3

model = Lar.extrudeSimplicial( VOID, repeat([1,-1],outer=10) )
Plasm.view(model)
model = Lar.extrudeSimplicial( model, repeat([1,-1],outer=10) )
Plasm.view(model)


# example 4

grid_2d = Lar.simplexGrid([3,3])
Plasm.view(grid_2d)
grid_3d = Lar.simplexGrid([2,3,4])
Plasm.view(grid_3d)
V,CV = Lar.simplexGrid([1,1,1])
Plasm.view(V,CV)

# example 5

SK2 = Lar.simplexFacets(CV)
Plasm.view(V, SK2)
SK1 = Lar.simplexFacets(SK2)
Plasm.view(V, SK1)

