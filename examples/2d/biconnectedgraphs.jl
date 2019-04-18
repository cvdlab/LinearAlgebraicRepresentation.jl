using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

n = 50 #3000 #1000 #20000
t = 0.75 #0.15 #0.4 #0.15

function randomlines(n=1000,t=0.2)
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
	return model
end

model = randomlines(n,t)
Sigma = Lar.spaceindex(model)
Plasm.view(model)

W,EW = Lar.fragmentlines(model)
Plasm.viewexploded(W,EW)(1.2,1.2,1.2)

V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)
hpcs = [ Plasm.lar2hpc(V,EVs[i]) for i=1:length(EVs) ]
Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

