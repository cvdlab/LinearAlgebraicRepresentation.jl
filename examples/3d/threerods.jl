using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using QHull

V,CV = Lar.rod()([4,1])
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


VV = [[k] for k=1:size(V,2)]
GL.VIEW( GL.numbering(.2)((V,[VV,EV,FV]), GL.COLORS[4]) );

cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');
V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)

W = convert(Lar.Points, V')
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)

GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
meshes = GL.GLExplode(V,CVs[1:end],8,4,6,99,1);
GL.VIEW( push!( meshes, GL.GLFrame) );
