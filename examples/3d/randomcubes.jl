using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL, LinearAlgebra
GL = ViewerGL

#include("")
store = [];
scaling = 1.5;
V,(VV,EV,FV,CV) = Lar.cuboid([0.25,0.25,0.25],true,[-0.25,-0.25,-0.25]);
mybox = (V,CV,FV,EV);

for k=1:10
	size = rand()*scaling
	scale = Lar.s(size,size,size)
	transl = Lar.t(rand(3)...)
	alpha = 2*pi*rand()
	rx = Lar.r(alpha,0,0); ry = Lar.r(0,alpha,0); rz = Lar.r(0,0,alpha)
	rot = rx * ry * rz
	str = Lar.Struct([ transl, scale, rot, mybox ])
	obj = Lar.struct2lar(str)
	vs = obj[1]
	diag = LinearAlgebra.norm(vs[:,8]-vs[:,1])
	if diag > 1/5
		push!(store, obj)
	end
end

str = Lar.Struct(store);
V,CV,FV,EV = Lar.struct2lar(str);
# V = Plasm.normalize3D(V) TODO:  solve MethodError bug

open("/tmp/lar.txt", "w") do f
	write(f, "V = $V\n\n")
	write(f, "CV = $CV\n\n")
	write(f, "FV = $FV\n\n")
	write(f, "EV = $EV\n\n")
	close(f)
end

GL.VIEW([ GL.GLPol(V,CV, GL.COLORS[2], 0.1) ]);

function testarrangement(V,CV,FV,EV)
		cop_EV = Lar.coboundary_0(EV::Lar.Cells);
		cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
		W = convert(Lar.Points, V');

		V, copEV, copFE, copCF = Lar.space_arrangement(
				W::Lar.Points, cop_EV::Lar.ChainOp, cop_FE::Lar.ChainOp);

		V = convert(Lar.Points, V');
		V,CVs,FVs,EVs = Lar.pols2tria(V, copEV, copFE, copCF) # whole assembly
		GL.VIEW(GL.GLExplode(V,FVs,1.1,1.1,1.1,99,1));
		GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
		GL.VIEW(GL.GLExplode(V,CVs,1,1,1,99,0.2));
end

testarrangement(V,CV,FV,EV);
