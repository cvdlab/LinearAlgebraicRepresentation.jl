using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

n = 20
t = 1/n
V = zeros(Float64,2,2*n)
EV = [zeros(Int64,2) for k=1:n]

for k=1:n
	v1 = rand(Float64,2)*t
	v2 = rand(Float64,2)*t
	vm = (v1+v2)/2
	V[:,k] = (v1-vm)*t + [rand(), rand()]
	V[:,2*k] = (v2-vm)*t + [rand(), rand()]
	EV[k] = [k,2*k]
end

Sigma = Lar.spaceindex((V,EV))

for k=1:length(Sigma)
	@show k," ",sort(Sigma[k])
end

Plasm.view(V,EV)


