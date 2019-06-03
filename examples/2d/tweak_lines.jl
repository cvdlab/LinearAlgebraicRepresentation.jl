using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

function extendlines(V::Lar.Points, EV::Lar.Cells, s=1.001)
	@assert size(V,1)==2 "extendlines() is for 2D models!"
	segments = [hcat(V[:,v1], V[:,v2]) for (v1,v2) in EV]
	n = length(EV)

	for segi in segments
		# renaming of segment vertices
		x1 = segi[1,1]; x2 = segi[1,2];
		y1 = segi[2,1]; y2 = segi[2,2];
		# middle segment point
		tx = sum(x1+x2)/2
		ty = sum(y1+y2)/2
		# composite translation and scaling transformation
		segi[1,1] = s*(x1-tx)+tx; segi[1,2] = s*(x2-tx)+tx;
		segi[2,1] = s*(y1-ty)+ty; segi[2,2] = s*(y2-ty)+ty;
	end
	V = zeros(Float64, 2, 2*length(EV))
	# new vertices for expanded line segments
	for (i,segi) in enumerate(segments)
		V[:,i] = segi[:,1]
		V[:,n+i] = segi[:,2]
	end
	EV = [[k,n+k] for k=1:length(segments)]
	return V, EV
end

n=30; t=2
V = zeros(Float64,2,2*n)
EV = [zeros(Int64,2) for k=1:n]

for k=1:n
	v1 = rand(Float64,2)
	v2 = rand(Float64,2)
	vm = (v1+v2)/2
	transl = rand(Float64,2)
	V[:,k] = (v1-vm)*t + transl
	V[:,n+k] = (v2-vm)*t + transl
	EV[k] = [k,n+k]
end

V = Plasm.normalize(V)
model = (V,EV)
Sigma = Lar.spaceindex(model)

Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))
Z,EZ = extendlines(V,EV)
Plasm.view(Plasm.numbering(.1)((Z,[[[k] for k=1:size(Z,2)], EZ])))

Plasm.view(V,EV)
model = V,EV;
model2 = Z,EZ;
W,EW = Lar.fragmentlines(model);
Z,EZ = Lar.fragmentlines(model2);

Plasm.viewexploded(W,EW)(1.2,1.2,1.2)
Plasm.viewexploded(Z,EZ)(1.2,1.2,1.2)
W,EW = Lar.fragmentlines(model);
Z,EZ = Lar.fragmentlines(model2);
U,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells));
T,ETs = Lar.biconnectedComponent((Z,EZ::Lar.Cells));
# 2-connected components (H & T)

W = convert(Lar.Points, U')
T = convert(Lar.Points, T')
cop_EV = Lar.coboundary_0(EW::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
cop_EZ = Lar.coboundary_0(EZ::Lar.Cells)
cop_EZ = convert(Lar.ChainOp, cop_EZ)
V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
Z, copEZ, copFE2 = Lar.Arrangement.planar_arrangement(T::Lar.Points, cop_EZ::Lar.ChainOp)

triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
triangulated_faces2 = Lar.triangulate2D(Z, [copEZ, copFE2])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
FZs = convert(Array{Lar.Cells}, triangulated_faces2)
W = convert(Lar.Points, V')
Z = convert(Lar.Points, Z')
Plasm.viewcolor(W::Lar.Points, FVs::Array{Lar.Cells})
Plasm.viewcolor(Z::Lar.Points, FZs::Array{Lar.Cells})

model = U,EVs;
model2 = T',ETs;
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
Plasm.view(Plasm.lar_exploded(model2)(1.2,1.2,1.2))
