using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

m1,M1 = minmax(rand(),rand())
m2,M2 = minmax(rand(),rand())
δ = 10.0^-2
v,(vv,ev,fv) = Lar.cuboid([1,1],true)

assembly = []
for k=1:10
	m1,M1 = minmax(rand(),rand())
	m2,M2 = minmax(rand(),rand())
	α = π * rand()
	β = 0.5 * rand()
	m = [(m1+M2)/2, (m2+M2)/2]
	V,(VV,EV,FV) = Lar.cuboid([M1+δ,M2+δ],true,[m1-δ,m2-δ])
	push!(assembly, Lar.Struct([Lar.t(m...)*Lar.r(α)*Lar.s(β,β)*Lar.t(-m...), (V,FV,EV)]))
	push!(assembly, (v,fv,ev))
end
assembly = Lar.struct2lar(Lar.Struct(assembly))
Plasm.view(assembly[1],assembly[3])
W,FW,EW = assembly

#for v=1:size(W,2)
#	W[:,v] = map(Lar.approxVal(4), W[:,v])
#end

V,(EV,FV),(copEV,copFE) = Lar.chaincomplex( W,EW )
Plasm.view(Plasm.numbering(.0625)((V,[[[k] for k=1:size(V,2)],EV,FV]))) 

c2 = sparseones
