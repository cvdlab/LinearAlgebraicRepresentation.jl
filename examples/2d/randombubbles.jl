using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

function randombubbles()
	store = []
	for k=1:50
		sc = rand()
		n = rand(3:48)
		scale = Lar.s(0.25*sc, 0.25*sc)
		transl = Lar.t(rand(2)...)
		str = Lar.Struct([ transl * scale, Lar.circle()() ])
		push!(store, str)
	end
	str = Lar.Struct(store);
	V,EV = Lar.struct2lar(str);
	V = GL.normalize2(V);
	V,FVs,EVs = Lar.arrange2D(V,EV)
end

# ////////////////////////////////////////////////////////////
# generation of 2D arrangement
GL.VIEW([ GL.GLLines(V,EV) ]);
V,FVs,EVs = randombubbles();

# native OpenGL visualization
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,1,1));
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));
GL.VIEW(GL.GLExplode(V,FVs,1.2,1.2,1.2,99,1));

GL.VIEW(GL.GLExplode(V,EVs,1.2,1.2,1.2,1,1));
