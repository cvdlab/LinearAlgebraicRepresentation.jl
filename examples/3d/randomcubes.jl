using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL, LinearAlgebra
GL = ViewerGL

#include("")
store = [];
scaling = .85;
V,(VV,EV,FV,CV) = Lar.cuboid([0.25,0.25,0.25],true,[-0.25,-0.25,-0.25]);
mybox = (V,CV,FV,EV);
for k=1:20
	size = rand()*scaling
	scale = Lar.s(size,size,size)
	transl = Lar.t(rand(3)...)
	alpha = 2*pi*rand()
	rx = Lar.r(alpha,0,0); ry = Lar.r(0,alpha,0); rz = Lar.r(0,0,alpha)
	rot = rx * ry * rz
	s = Lar.Struct([ transl, scale, rot, mybox ])
	obj = Lar.struct2lar(s)
	vs = obj[1]
	diag = LinearAlgebra.norm(vs[:,8]-vs[:,1])
	if diag > 1/10
		push!(store, obj)
	end
end

s = Lar.Struct(store);
V,CV,FV,EV = Lar.struct2lar(s);
# V = Plasm.normalize3D(V) TODO:  solve MethodError bug

GL.VIEW([ GL.GLPol(V,CV, GL.COLORS[1]) ]);


cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement(
	W::Lar.Points, cop_EV::Lar.ChainOp, cop_FE::Lar.ChainOp);


triangulated_faces = Lar.triangulate2D(V, [copEV, copFE]);
FVs = convert(Array{Lar.Cells}, triangulated_faces);
V = convert(Lar.Points, V');
GL.VIEW( GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1) );


EVs = Lar.FV2EVs(copEV, copFE); # polygonal face fragments
GL.VIEW( GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1) );
