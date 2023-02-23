using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function twosquares()
    V, (VV, EV, FV) = Lar.cuboid([1,1], true)
    assembly = Lar.Struct([ (V, EV, FV), Lar.t(.25,.25), (V, EV, FV) ])
    V,EW,FW = Lar.struct2lar(assembly)
    W = convert(Lar.Points, V')
    cop_EV = Lar.lar2cop(EW)
    V, copEV, copFE = Lar.Arrangement.planar_arrangement(W, cop_EV);
    triangulated_faces = Lar.triangulate2D(V, [copEV, copFE]);
    FVs = convert(Array{Lar.Cells}, triangulated_faces);
    W = convert(Lar.Points, V');
    GL.VIEW(GL.GLExplode(W,FVs,1.2,1.2,1.2,99,1));
    GL.VIEW(GL.GLExplode(W,FVs,1,1,1,99,1));
end

twosquares()

