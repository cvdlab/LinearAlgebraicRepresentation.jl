using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using QHull

V,CV = Lar.rod()()
ch = QHull.chull(convert(Lar.Points,V'))
FV = ch.simplices
EV = Lar.simplexFacets(FV)
GL.VIEW([ GL.GLGrid(V,EV) ]);

cyl = Lar.Struct([ Lar.t(0,0,-1.5), (V,FV,EV) ])
tris = Lar.Struct([ cyl,
                    Lar.Struct([Lar.r(pi/2,0,0), cyl]),
                    Lar.Struct([Lar.r(0,pi/2,0), cyl])
                    ])
V,FV,EV = Lar.struct2lar(tris)
GL.VIEW([ GL.GLGrid(V,FV) ]);
