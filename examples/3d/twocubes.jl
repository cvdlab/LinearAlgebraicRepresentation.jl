using LinearAlgebraicRepresentation
using Plasm, SparseArrays
Lar = LinearAlgebraicRepresentation
using Revise

function twocubes()
    #V,(VV,EV,FV,CV) = Lar.cuboid([0.5,0.5,0.5],true,[-0.5,-0.5,-0.5])
    V,(VV,EV,FV,CV) = Lar.cuboidGrid([4,4,4],true)
    mybox = (V,CV,FV,EV)
    #mybox = Lar.Struct([ Lar.t(0,0,1.), mybox ])

    twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), mybox ])
    #twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/3,0,0), Lar.r(0,0,pi/6), mybox ])
    V,CV,FV,EV = Lar.struct2lar(twocubes)
    Plasm.view(V,CV)

    cop_EV = Lar.coboundary_0(EV::Lar.Cells);
    cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
    W = convert(Lar.Points, V');

    V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W::Lar.Points, cop_EV::Lar.ChainOp, cop_FE::Lar.ChainOp)

    EV = Lar.cop2lar(copEV)
    FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
    FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]
    FV = convert(Lar.Cells, FV)
    W = convert(Lar.Points, V')

    #Plasm.view(Plasm.numbering(0.25)((W,[[[k] for k=1:size(W,2)],EV,FV])))
    Plasm.view(W,EV)

    triangulated_faces = Lar.triangulate(V, [copEV, copFE])
    FVs = convert(Array{Lar.Cells}, triangulated_faces)
    V = convert(Lar.Points, V')
    Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

    EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
    model = V,EVs
    Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
end

twocubes()
