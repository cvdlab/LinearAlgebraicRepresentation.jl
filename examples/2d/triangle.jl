using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
V,(VV,EV,FV) = Lar.simplex(2, true)
triangle = (V,EV,FV)
model = Lar.Struct([ triangle, Lar.t(.15,.15), Lar.s(.5,.5), triangle ])
V,EV,FV = Lar.struct2lar(model)
Plasm.view(Plasm.numbering(0.5)((V,[[[k] for k=1:size(V,2)], EV])))
