using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL


function randlines(n=300, t=0.4)
	#n = 100 #1000 #1000 #20000
	#t = 0.4 #0.15 #0.4 #0.15
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
	Sigma = Lar.spaceindex(model)

	model = V,EV;
	W,EW = Lar.fragmentlines(model);
	U,EVs = Lar.biconnectedComponent((W,EW::Lar.Cells));
	EV = convert(Lar.Cells, cat(EVs))
	V,FVs,EVs = Lar.arrange2D(U,EV)
end


# ////////////////////////////////////////////////////////////

# generation of 2D arrangement
V,FVs,EVs = randlines()

# native OpenGL visualization
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,3));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99));
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99));
