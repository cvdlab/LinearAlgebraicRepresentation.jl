using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm


function cuboids(n,scale=1.)
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		dims = rand(Float64, 2)
		V,(_,EV,_) = Lar.cuboid(corner,true,corner+dims)
		center = (corner + corner+dims)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end



V,EV = cuboids(50, .15)
V = Plasm.normalize(V,flag=true)
model2d = V,EV

Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))

Sigma =  Lar.spaceindex(model2d);
println([(k,Sigma[k]) for k=1:length(Sigma)])

Plasm.view(Plasm.numbering(.1)((V,[[[k] for k=1:size(V,2)], EV])))



