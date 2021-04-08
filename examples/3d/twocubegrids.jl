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
    GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1], 0.5) ]);

    cop_EV = Lar.coboundary_0(EV::Lar.Cells);
    cop_EW = convert(Lar.ChainOp, cop_EV);
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    W = convert(Lar.Points, V');

    V, copEV, copFE, copCF = Lar.space_arrangement( W, cop_EW, cop_FE)

    EV = Lar.cop2lar(copEV)
    FE = Lar.cop2lar(copFE)
    CF = Lar.cop2lar(copCF)
    FV = [union([EV[e] for e in fe]...) for fe in FE]
    CV = [union([FV[f] for f in cf]...) for cf in CF]
    W = convert(Lar.Points, V')
    WW = [[k] for k=1:size(W,2)]

    GL.VIEW(GL.numbering(.5)((W,[WW,EV,FV,CV])));

    V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF)
    GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99,1));
    GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,99,1));
    GL.VIEW(GL.GLExplode(V,CVs[2:end],5,5,5,99,0.5));

    return V, copEV, copFE, copCF
end

twocubegrids(3,3,3);
