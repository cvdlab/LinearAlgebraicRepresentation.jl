using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm


function randomcuboids(n,scale=1.)
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



V,EV = randomcuboids(100, .2)
V = Plasm.normalize(V,flag=true)
model2d = V,EV

Sigma =  Lar.spaceindex(model2d);
for k=1:length(Sigma) println(k,Sigma[k]) end

Plasm.view(Plasm.numbering(.15)((V,[[[k] for k=1:size(V,2)], EV])))
Plasm.viewexploded(V,EV)(1.2,1.2,1.2) 	# no numerical errors

W,EW = Lar.fragmentlines((V,EV)) 
Plasm.viewexploded(W,EW)(1.2,1.2,1.2)	
Plasm.view(Plasm.numbering(.015)((W,[[[k] for k=1:size(W,2)], EW])))

V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)
hpcs = [ Plasm.lar2hpc(V,EVs[i]) for i=1:length(EVs) ]
Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

FVs = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewlarcolor(V::Lar.Points, FVs::Array{Lar.Cells})

