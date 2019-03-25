using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm


function cuboids(n,scale=1.)
	assembly = []
	for k=1:n
		corner = rand(Float64, 2)
		sizes = rand(Float64, 2)
		V,(_,EV,_) = Lar.cuboid(corner,true,corner+sizes)
		center = (corner + corner+sizes)/2
		angle = rand(Float64)*2*pi
		obj = Lar.Struct([ Lar.t(center...), Lar.r(angle), 
				Lar.s(scale,scale), Lar.t(-center...), (V,EV) ])
		push!(assembly, obj)
	end
	Lar.struct2lar(Lar.Struct(assembly))
end



V,EV = cuboids(20, .35)
V = Plasm.normalize(V,flag=true)
model2d = V,EV

Sigma =  Lar.spaceindex(model2d);
for k=1:length(Sigma) println(k,Sigma[k]) end

Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV])))




