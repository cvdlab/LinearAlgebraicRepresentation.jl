using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function randombubbles()
	store = []
	mycircle(r,n) = Lar.circle(r)([n])
	for k=1:50
		r = rand()
		n = 8
		while true
			n = abs(rand(Int8)+1)
			if n>2 break end
		end
		scale = Lar.s(0.25,0.25)
		transl = Lar.t(rand(2)...)
		s = Lar.Struct([ transl, scale, mycircle(r,n) ])
		push!(store, Lar.struct2lar(s))
	end
	s = Lar.Struct(store)
	V,EV = Lar.struct2lar(s)
	V = GL.normalize2(V)
	GL.VIEW([ GL.GLLines(V,EV) ])

	W,EW = Lar.fragmentlines((V,EV))
	V,FVs,EVs = Lar.arrange2D(W,EW)
end

# ////////////////////////////////////////////////////////////
# generation of 2D arrangement
V,FVs,EVs = randombubbles()

# native OpenGL visualization
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,1,1));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,3,1));
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99,1));

GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,1,1));
GL.VIEW(GL.GLExplode(V,EVs,1,1,1,1,1));
