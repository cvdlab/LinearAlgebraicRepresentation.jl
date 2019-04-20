using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation

V,(VV,EV,FV,CV) = Lar.cuboid([0.5,0.5,0.5],true,[-0.5,-0.5,-0.5])
mybox = (V,CV,FV,EV)

twocubes = Lar.Struct([ mybox, Lar.t(0.3,0.4,0.5), Lar.r(pi/3,0,0), Lar.r(0,0,pi/6), mybox ])
V,CV,FV,EV = Lar.struct2lar(twocubes)
Plasm.view(V,CV)

cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');
V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement(
	W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)

