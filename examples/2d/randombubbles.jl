using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

function randombubbles()
	store = []
	mycircle(r,n) = Lar.circle(r)(n)
	for k=1:30
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
	V = Plasm.normalize(V)
	#Plasm.view(V,EV)
	GL.VIEW([ GL.GLLines(V,EV) ])
	V,FVs = Lar.arrange2D(V,EV)
end

# ////////////////////////////////////////////////////////////

# generation of 2D arrangement
V,FVs = randombubbles()

# native OpenGL visualization
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,3));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99));
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99));
