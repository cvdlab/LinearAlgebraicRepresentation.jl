using LinearAlgebraicRepresentation
using SparseArrays
Lar = LinearAlgebraicRepresentation
L = Lar
using ViewerGL
GL = ViewerGL

function twocubegrids(n,m,p)
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
    mybox = (V,CV,FV,EV)

    twocubs = Lar.Struct([mybox, L.t(.3,.4,.5), L.r(pi/5,0,0), L.r(0,0,pi/12), mybox])

    V,CV,FV,EV = Lar.struct2lar(twocubs)
    GL.VIEW([ GL.GLGrid(V,FV, GL.Point4d(1,1,1,0.2)) ]);

    cop_EV = Lar.coboundary_0(EV::Lar.Cells);
    cop_EW = convert(Lar.ChainOp, cop_EV);
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    W = convert(Lar.Points, V');

    V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EW, cop_FE)

    EV = Lar.cop2lar(copEV)
    FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
    FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]
    FV = convert(Lar.Cells, FV)
    W = convert(Lar.Points, V')
    WW = [[k] for k=1:size(W,2)]

    GL.VIEW(GL.numbering(.2)((W,[WW,EV])));

    triangulated_faces = Lar.triangulate(V, [copEV, copFE])
    FVs = convert(Array{Lar.Cells}, triangulated_faces)
    V = convert(Lar.Points, V')
    GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));
    EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
    GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));

    return V, copEV, copFE, copCF
end

twocubegrids(5,5,5);
