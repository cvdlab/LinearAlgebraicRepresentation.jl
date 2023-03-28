# polyhedron visualization and integration test
# ---------------------------------------------
using ViewerGL
using LinearAlgebraicRepresentation
GL = ViewerGL
Lar = LinearAlgebraicRepresentation

# your vertices
W = [-1.0 1.0 1.0
    -1.0 -1.0 1.0
    1.0 0.0 -1.0
    1.0 1.0 1.0
    0.0 1.0 -1.0]
    
# Transform in our vertex repr (transpose)
W = convert(Matrix,W')

# added two more triangle faces (your polyhedron is not linear)
# it has a face with four non-coplanar points (my last two triangles)
TW = [[1,2,4],[2,3,4],[3,5,4],[1,4,5],[2,1,5],[2,5,3]]

# LAR representation of 8 unit cubes in positive octant
V,(VV,EV,FV,CV) = Lar.cuboidGrid([2,2,2],true)

# translate your polyhedral model to have the conic vertex in [2,2,2]'
W,TW = Lar.apply(Lar.t(1,1,1))((W,TW))

# you may use the mouse to rotate/scale
GL.VIEW([
   GL.GLGrid(V,EV,GL.COLORS[1],1.)
   GL.GLGrid(W,TW,GL.COLORS[6],0.5)
   GL.GLAxis(GL.Point3d(1,1,1),GL.Point3d(2,2,2))
]);

# computation of 4x4 affine matrix of inertia and volume/surface using boundary triangulation, according to:

# C. Cattani and A. Paoluzzi: Boundary integration over linear polyhedra. Computer-Aided Design 22(2): 130-135 (1990) (doi:10.1016/0010-4485(90)90007-Y) 

pol = (W,TW)
vol = Lar.volume(pol)