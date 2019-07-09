using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

n,m,p = 2,2,2
V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
cube = V,FV,EV

threecubes = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube ])

V,FV,EV = Lar.struct2lar(threecubes)
GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame ])

cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE)

#triangulated_faces = Lar.triangulate(V, [copEV, copFE])
W = convert(Lar.Points, V')
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)

GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
meshes = GL.GLExplode(V,CVs[1:end],8,4,6,99,0.5);
GL.VIEW( push!( meshes, GL.GLFrame) );
