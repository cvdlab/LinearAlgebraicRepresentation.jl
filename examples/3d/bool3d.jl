using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

# 3D Boolean example generation
#-------------------------------------------------------------------------------
n,m,p = 1,1,1
V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
cube = V,FV,EV

threecubes = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube ])

V,FV,EV = Lar.struct2lar(threecubes)
GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame ]);

boolmatrix = Lar.booleanops3d(threecubes)
Matrix(boolmatrix)
#three-chains = [ for k = 1:3]
