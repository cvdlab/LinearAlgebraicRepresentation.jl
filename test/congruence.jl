using LinearAlgebraicRepresentation
using ViewerGL, LinearAlgebra, SparseArrays, DataStructures, NearestNeighbors
L = Lar = LinearAlgebraicRepresentation
GL = ViewerGL

# initialization
n = 3
v,(_,ev,fv) = L.cuboidGrid([n,n],true)
v = [v; zeros(1,(n+1)*(n+1))]
g = (v,ev,fv)

# range di numeri casuali piccoli 0e10-2  x  1e10-2
interval =  0e10-2 : 1e10-2

# matrice 4x4 di numeri casuali piccoli (dipende dagli zeri)
A(interval) = (1_000_000 ./ rand(interval,4,4))

# Esempio di complesso 2D immerso in 3D
G1 = L.Struct([  g, L.r(pi/2,0,0) + A(interval),  g , L.r(0,pi/2,0), g , L.r(0,0,pi/2), g , 
		L.t(2,0,0), L.r(pi/2,0,0),  g , L.t(0,0,2), L.r(0,pi/2,0), g ]);
V1,EV1,FV1 = Lar.struct2lar(G1);
VV1 = [[k] for  k=1:size(V1,2)];
GL.VIEW(push!(GL.numbering(.5)((V1,[VV1,EV1,FV1]), GL.COLORS[1], 0.1),GL.GLFrame2));

# Esempio di complesso 2D immerso in 3D
G2 = L.Struct([  g, L.r(pi/2,0,0) + A(interval)./10, g , L.r(0,pi/2,0), g , L.r(0,0,pi/2), g , 
		L.t(2,0,0), L.r(pi/2,0,0),  g , L.t(0,0,2), L.r(0,pi/2,0), g ]);
V2,EV2,FV2 = Lar.struct2lar(G2);
VV2 = [[k] for  k=1:size(V2,2)];
GL.VIEW(push!(GL.numbering(.5)((V2,[VV2,EV2,FV2]), GL.COLORS[1], 0.1),GL.GLFrame2));

# complesso unione combinatoria di G1 e G2
G = L.Struct([ G1, G2 ])
V,EV,FV = Lar.struct2lar(G)
VV = [[k] for  k=1:size(V,2)]
GL.VIEW(push!(GL.numbering(.5)((V,[VV,EV,FV]), GL.COLORS[1], 0.1),GL.GLFrame2));
Delta_0 = Lar.coboundary_0(EV::Lar.Cells)
Delta_1 = Lar.coboundary_1(FV::Lar.Cells, EV::Lar.Cells)

W,EW,FE,FW = Lar.chaincongruence(V, Delta_0, Delta_1, 1e-1);

WW = [[k] for  k=1:size(W,2)]
GL.VIEW(push!(GL.numbering(.5)((Lar.Points(W),Lar.Cells[WW,EW,FW]), GL.COLORS[1], 0.1),GL.GLFrame2));
