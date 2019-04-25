using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

function randomlines(n=1000, t=0.3)
	#n = 1000 #1000 #1000 #20000
	#t = 0.30 #0.15 #0.4 #0.15
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
	#Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

	Plasm.view(V,EV)
	model = V,EV;
	W,EW = Lar.fragmentlines(model);

	Plasm.viewexploded(W,EW)(1.2,1.2,1.2)
	W,EW = Lar.fragmentlines(model);
	V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells));
	# 2-connected components (H & T)

	model = V,EVs;
	Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
end

randomlines()
