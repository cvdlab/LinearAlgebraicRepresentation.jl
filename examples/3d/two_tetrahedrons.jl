using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

V, (VV,EV,FV,CV) = Lar.simplex(3, true)
tetra = V, EV,FV,CV
twotetra = Lar.Struct([ tetra, Lar.t(0.25,0.25,0.25), tetra ])
V,EV,FV,CV = Lar.struct2lar(twotetra)
GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1],.2) ]);

GL.VIEW(GL.numbering(.5)((V,[[[v] for v=1:size(V,2)],EV,FV,CV]) ));

cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

V, copEV, copFE, copCF = Lar.space_arrangement( W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)

# copCF = OK !!

V,CVs,FVs,EVs = Lar.pols2tria(V, copEV, copFE, copCF); # whole assembly
GL.VIEW(GL.GLExplode(V,FVs,1.1,1.1,1.1,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
sort!( CVs; by=x->length(x) );
GL.VIEW(GL.GLExplode(V,CVs[1:end-1]
    ,1.,1.,1.,99,0.5));
GL.VIEW(GL.GLExplode(V,[CVs[end]],4.,4.,4.,99,0.5));

#EV = Lar.cop2lar(copEV)
#CV = [collect(Set(Lar.cat([FV[f]  for  f in cell]))) for cell in cop2lar(copCF)]
#FE = [findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
#FV = [union(EV[e] for e in FE[f]) for f=1:length(FE)]
#FV = convert(Lar.Cells, FV)
#W = convert(Lar.Points, V') # V'
#WW = [[k] for k=1:size(W,2)]
#
#
#triangulated_faces = Lar.triangulate(V, [copEV, copFE])
#FVs = convert(Array{Lar.Cells}, triangulated_faces)
#V = convert(Lar.Points, V')
#GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));
#
#EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
#GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
