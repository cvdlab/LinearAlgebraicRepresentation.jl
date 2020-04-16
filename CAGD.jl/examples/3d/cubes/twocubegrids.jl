todisplay = VERSION <= VersionNumber("1.2") ? true : false

using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
if todisplay
    using ViewerGL
    GL = ViewerGL
	include("../views.jl")
end

function twocubegrids(n=1, m=1, p=1, atol = 1e-6)
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
    mybox = (V,CV,FV,EV)

    twocubs = Lar.Struct([mybox, Lar.t(.3,.4,.5), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), mybox])

    V,CV,FV,EV = Lar.struct2lar(twocubs)
    if todisplay  GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1]) ])  end

    cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    model = CAGD.Model(V)
    CAGD.addModelCells!(model, 1, cop_EV)
    CAGD.addModelCells!(model, 2, cop_FE)

    split_model = CAGD.pairwise_decomposition(model)
    congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)
    if todisplay  displayModel(congr_model)  end
    gift_model, bicon_comps = CAGD.tgw(congr_model, 3, atol = atol)
    if todisplay  viewExplode(gift_model)  end
end