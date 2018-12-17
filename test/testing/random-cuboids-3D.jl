using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

m1,M1 = minmax(rand(),rand())
m2,M2 = minmax(rand(),rand())
m3,M3 = minmax(rand(),rand())
δ = 10.0^-2

assembly = []
for k=1:10
	m1,M1 = minmax(rand(),rand())
	m2,M2 = minmax(rand(),rand())
	m3,M3 = minmax(rand(),rand())
	α = π * rand()
	m = [(m1+M2)/2, (m2+M2)/2, (m3+M3)/2]
	V,(VV,EV,FV,CV) = Lar.cuboid([M1+δ,M2+δ,M3+δ],true,[m1-δ,m2-δ,m3-δ])
	push!(assembly, #Lar.Struct([
		Lar.t(m...)*Lar.r(0,0,α)*Lar.t(-m...), 
		(V,FV,EV)
		#])
		)
end
assembly = Lar.struct2lar(Lar.Struct(assembly))
Plasm.view(assembly[1:2])
W,FW,EW = assembly

for v=1:size(W,2)
	W[:,v] = map(Lar.approxVal(4), W[:,v])
end

V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW );
@show CV
Plasm.view(Plasm.numbering(.0625)((V,[[[k] for k=1:size(V,2)],EV,FV]))) 

cc = [copEV,copFE,copCF]
cc = convert(Array{SparseMatrixCSC{Int,Int64},1}, cc)
Lar.lar2obj(V'::Lar.Points, cc)
