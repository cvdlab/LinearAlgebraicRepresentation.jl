using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

n = 200 #3000 #1000 #20000
t = 0.35 #0.15 #0.4 #0.15

function randlines(n=1000,t=0.2)
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
	V = GL.normalize2(V)
	model = (V,EV)
	return model
end

model = randlines(n,t)
Sigma = Lar.spaceindex(model)
GL.VIEW([ GL.GLLines(model...)] )

W,EW = Lar.fraglines(1.5,1.5,1.5)(model)
GL.VIEW([ GL.GLLines(W,EW,GL.COLORS[1])] );

W,EW = Lar.fragmentlines(model)
V,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells)) # 2-connected components (H & T)

graphs = [ GL.GLLines(V,EVs[k],GL.COLORS[k%12]) for k=1:length(EVs) ]
GL.VIEW(graphs);
