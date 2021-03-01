#=
1. Construct N parallel cuboids of size 10000 × 1 × 100 spaced one unit apart in y -direction as object W . 
2. Construct N parallel cuboids of size 1 × 10000 × 100 spaced one unit apart in x -direction as object W ′ . 
3. Align W and W ′ at their lower front left corner.
4. Move W ′ along the z-axis for fifty units.
5. Measure time for W ∪ W′.
=#

using ViewerGL; GL = ViewerGL
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation;

n = 30
V1,(_,EV1,FV1,_) = Lar.cuboid([1,2n,10],true);
cuboid1 = (V1,FV1,EV1);
V2,(_,EV2,FV2,_) = Lar.cuboid([2n,1,10],true);
cuboid2 = (V2,FV2,EV2);

# @time let 
	group1 = Lar.struct2lar(Lar.Struct(repeat([cuboid1, Lar.t(2,0,0)],outer=n)));
	group2 = Lar.struct2lar(Lar.Struct(repeat([cuboid2, Lar.t(0,2,0)],outer=n)));
	assembly = Lar.Struct([group1, Lar.t(0,0,5), group2]);
	
	V,FV,EV = Lar.struct2lar(assembly);
	cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
	cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells); ## TODO: debug
	W = convert(Lar.Points, V');
	V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W, cop_EV, cop_FE);
	W = convert(Lar.Points, V');
	V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF);
# end;
	
GL.VIEW( GL.GLExplode(V,CVs[2:end],5,5,5,99,0.5) );
GL.VIEW( GL.GLExplode(V,CVs[1:1],5,5,5,99,0.5) );
