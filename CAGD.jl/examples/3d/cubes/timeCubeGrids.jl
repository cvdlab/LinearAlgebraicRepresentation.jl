using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using BenchmarkTools
using CAGD

function generateCubeGrids(n=1, m=1, p=1)
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
    mybox = (V,CV,FV,EV)

    twocubs = Lar.Struct([mybox, Lar.t(.3,.4,.5), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), mybox])

    V,CV,FV,EV = Lar.struct2lar(twocubs)

    cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    model = CAGD.Model(V)
    CAGD.addModelCells!(model, 1, cop_EV)
    CAGD.addModelCells!(model, 2, cop_FE)

    return model
end

m = generateCubeGrids(1, 1, 1)

@btime CAGD.spatial_arrangement(m)