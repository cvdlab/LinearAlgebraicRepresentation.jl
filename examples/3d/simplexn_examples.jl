using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL


VOID = [[]], [[1]]


# example 1

V = [[0,0] [1,0] [2,0] [0,1] [1,1] [2,1] [0,2] [1,2] [2,2]];
FV = [[1,2,4],[2,3,5],[3,5,6],[4,5,7],[5,7,8],[6,8,9]];
pattern = repeat([1,2,-3],outer=4);
model = (V,FV);
W,FW = Lar.extrudeSimplicial(model, pattern);
GL.VIEW([ GL.GLGrid(W,FW, GL.COLORS[7], 0.5) ]);

# example 2

model = Lar.extrudeSimplicial( VOID, ones(10) )
GL.VIEW([ GL.GLGrid(model..., GL.COLORS[7], 0.5) ]); #TODO: check

model = Lar.extrudeSimplicial( model, ones(10) )
GL.VIEW([ GL.GLLines(model..., GL.COLORS[7]),
        GL.GLAxis(GL.Point3d(0,0,0),GL.Point3d(1,1,1)) ]);

model = Lar.extrudeSimplicial( model, ones(10) )
GL.VIEW([ GL.GLGrid(model..., GL.COLORS[7], 0.5) ]);

# example 3

model = Lar.extrudeSimplicial( VOID, repeat([1,-1],outer=10) )
GL.VIEW([ GL.GLGrid(model..., GL.COLORS[7], 0.5) ]); #TODO: check

model = Lar.extrudeSimplicial( model, repeat([1,-1],outer=10) )
GL.VIEW([ GL.GLGrid(model..., GL.COLORS[7],0.7),
        GL.GLAxis(GL.Point3d(0,0,0),GL.Point3d(1,1,1)) ]);


# example 4

grid_2d = Lar.simplexGrid([3,3])
GL.VIEW([ GL.GLGrid(grid_2d..., GL.COLORS[7],0.7),
        GL.GLLines(grid_2d...),    # TODO:  check lines
        GL.GLAxis(GL.Point3d(0,0,0),GL.Point3d(1,1,1)) ]);

grid_3d = Lar.simplexGrid([2,3,4])
GL.VIEW([ GL.GLGrid(grid_3d..., GL.COLORS[7]),
        GL.GLAxis(GL.Point3d(0,0,0),GL.Point3d(1,1,1)) ]);  # TODO  check
