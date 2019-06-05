using LinearAlgebraicRepresentation
using Plasm, LinearAlgebra
Lar = LinearAlgebraicRepresentation
using Revise, SparseArrays

store = []
scaling = 1.75
V,(VV,EV,FV,CV) = Lar.cuboid([0.25,0.25,0.25],true,[-0.25,-0.25,-0.25])
mybox = (V,CV,FV,EV)
for k=1:5
	size = rand()*scaling
	scale = Lar.s(size,size,size)
	transl = Lar.t(rand(3)...)
	alpha = 2*pi*rand()
	rx = Lar.r(alpha,0,0); ry = Lar.r(0,alpha,0); rz = Lar.r(0,0,alpha)
	rot = rx * ry * rz
	s = Lar.Struct([ transl, scale, rot, mybox ])
	push!(store, Lar.struct2lar(s))
end

s = Lar.Struct(store)
V,CV,FV,EV = Lar.struct2lar(s)
# V = Plasm.normalize3D(V) TODO:  solve MethodError bug
Plasm.view(V,CV)

cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

@show W;
@show findnz(cop_EW);
@show findnz(cop_FE);

V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement(
	W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)

cc = [copEV, copFE, copCF]
output = Lar.lar2obj(V::Lar.Points, cc::Lar.ChainComplex)

f = open("test/out3d.obj", "w")
print(f, output)
close(f)

V,EVs,FVs = Lar.obj2lar("test/out3d.obj")

for k=1:length(FVs)
	Plasm.view(V, FVs[k])
end
