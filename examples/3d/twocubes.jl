using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL, SparseArrays
GL = ViewerGL

function twocubes()
    V, (VV, EV, FV, CV) = Lar.cuboid([1,1,1], true)
    assembly = Lar.Struct([ (V, EV, FV), Lar.t(0.5,0.5,.5), (V, EV, FV) ])
    V,EW,FW = Lar.struct2lar(assembly)
    W = convert(Lar.Points, V')
        
    cop_EW = Lar.lar2cop(EW)
    cop_FE = Lar.coboundary_1(V, FW::Lar.Cells, EW::Lar.Cells)
    V, copEV, copFE, copCF = Lar.space_arrangement(
        W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp);
	V,CVs,FVs,EVs = Lar.pols2tria(V', copEV, copFE, copCF); # whole assembly
	GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,1));
	GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
	GL.VIEW(GL.GLExplode(V,CVs[2:end],1.5,1.5,1.5,99,1));
end

twocubes()
