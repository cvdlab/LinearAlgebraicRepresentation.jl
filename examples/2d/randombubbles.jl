using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation

function randombubbles()
	store = []
	mycircle(r,n) = Lar.circle(r)(n)
	for k=1:30
		r = rand()
		global n = 8
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
	Plasm.view(V,EV)

	W = convert(Lar.Points, V')
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)
	V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

	triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	V = convert(Lar.Points, V')
	Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

	EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
	model = V,EVs
	Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
end

randombubbles()