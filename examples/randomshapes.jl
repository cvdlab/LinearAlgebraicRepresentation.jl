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



V,EV = cuboids(40, .35)
V = Plasm.normalize(V,flag=true)
model2d = V,EV

Sigma =  Lar.spaceindex(model2d);
for k=1:length(Sigma) println(k,Sigma[k]) end

Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV])))
Plasm.viewexploded(V,EV)(1.2,1.2,1.2) 	# no numerical errors

W,EW = Lar.fragmentlines((V,EV)) # introduced false edges (measure ≈ 0)
#W,EW = Lar.simplifyCells(W,EW)
Plasm.viewexploded(W,EW)(1.2,1.2,1.2)	
Plasm.view(Plasm.numbering(.0035)((W,[[[k] for k=1:size(W,2)], EW])))



