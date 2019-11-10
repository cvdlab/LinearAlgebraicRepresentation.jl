using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using QHull, Plasm

V,CV = Lar.rod()([4,1])
ch = QHull.chull(convert(Lar.Points,V'))
FV = ch.simplices
EV = Lar.simplexFacets(FV)
GL.VIEW([ GL.GLGrid(V,EV,GL.COLORS[1],1) ]);

cyl = Lar.Struct([ Lar.t(0,0,-1.5), (V,FV,EV) ])
cyl = Lar.struct2lar(cyl)

V,FV = ([1.0 1.0 -0.5 -0.5 -0.5 -0.5; 0.0 0.0 0.86603 0.86603 -0.86603 -0.86603; 0.0 3.0 0.0 3.0 0.0 3.0],
Array{Int64,1}[[3, 4, 1, 2], [3, 4, 5, 6], [5, 6, 1, 2], [3, 1, 5], [2, 4, 6]])
EV = [[1,3],[3,5],[5,1], [2,4],[4,6],[6,2], [1,2],[3,4],[5,6]]
cyl = (V,FV,EV)

tris = Lar.Struct([ cyl,
                    Lar.Struct([Lar.r(pi/2,0,0), cyl ]) ,
                    Lar.Struct([Lar.r(0,pi/2,0), cyl ]) ,
                    Lar.Struct([Lar.r(0,pi,0), cyl ]) ,
                    ])
V,FV,EV = Lar.struct2lar(tris)
Plasm.view(V,EV) #GL.VIEW([ GL.GLGrid(V,FV) ]);

#VV = [[k] for k=1:size(V,2)]
#GL.VIEW( GL.numbering(.2)((V,[VV,EV,FV]), GL.COLORS[4]) );

cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
#cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
cop_FE = Lar.coboundary_1(FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');
V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)

W = convert(Lar.Points, V')
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)

GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,99,1));
meshes = GL.GLExplode(V,CVs[1:end],4,4,4,99,0.5);
GL.VIEW( push!( meshes, GL.GLFrame) );
