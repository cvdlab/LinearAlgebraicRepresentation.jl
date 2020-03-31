display = VERSION < VersionNumber("1.2") ? true : false

using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
CAGD = CAGD
if display
    using ViewerGL
    GL = ViewerGL
end


function viewExplode(model)
    V,CVs,FVs,EVs = Lar.pols2tria(model.G, model.T[1], model.T[2], model.T[3])

    GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
    GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
    GL.VIEW(GL.GLExplode(V,CVs[2:end],5,5,5,99,0.5));
end

function displayModel(model)
    V = model.G
    EV = Lar.cop2lar(model.T[1])
    FE = Lar.cop2lar(model.T[2])
    CF = Lar.cop2lar(model.T[3])
    FV = Lar.cop2lar(map(x -> Int8(x/2), abs.(model.T[2]) * abs.(model.T[1])))
    CV = []
    triangulated_faces = Lar.triangulate(convert(Lar.Points, V'), [model.T[1], model.T[2]])
    FVs = convert(Array{Lar.Cells}, triangulated_faces)
    EVs = Lar.FV2EVs(model.T[1], model.T[2])

    GL.VIEW([ GL.GLPol(V,EV, GL.COLORS[1]) ]);
    GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));
    GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
end

function twocubegrids(n=1, m=1, p=1)
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
    mybox = (V,CV,FV,EV)

    twocubs = Lar.Struct([mybox, Lar.t(.3,.4,.5), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), mybox])

    V,CV,FV,EV = Lar.struct2lar(twocubs)
    if display  GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1]) ])  end

    cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells))
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    model = CAGD.Model(V)
    CAGD.addModelCells!(model, 1, cop_EV)
    CAGD.addModelCells!(model, 2, cop_FE)

    split_model = CAGD.pairwise_decomposition(model)
    congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)
    if display  displayModel(congr_model)  end
    gift_model, bicon_comps = CAGD.tgw(congr_model, 3)
    if display  viewExplode(gift_model)  end
end