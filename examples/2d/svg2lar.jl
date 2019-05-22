using Plasm,SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
import Base.show

function show(filename)
	V, EV = Plasm.svg2lar(filename)
	Plasm.view(V,EV)
	return V, EV
end

#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/new.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/curved.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/twopaths.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/paths.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/boundarytest2.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/tile.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/interior.svg")
#show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/holes.svg")

V,EV = show("/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/Lar.svg")

# subdivision of input edges
W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV, copFE = Lar.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

# compute containment graph of components
bicon_comps = Lar.Arrangement.biconnected_components(copEV)

# visualization of component graphs
EW = Lar.cop2lar(copEV)
W = convert(Lar.Points, V')
hpcs = [ Plasm.lar2hpc(W,[EW[e] for e in comp]) for comp in bicon_comps ]
Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

# final solid visualizations
FE = [SparseArrays.findnz(copFE[k,:])[1] for k=1:size(copFE,1)]
FV = [collect(Set(cat([EW[e] for e in FE[f]]))) for f=1:length(FE)]

triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
V = convert(Lar.Points, V')
FVs = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

EVs = Lar.FV2EVs(copEV, copFE) # polygonal face boundaries
EVs = convert(Array{Array{Array{Int64,1},1},1}, EVs)
Plasm.viewcolor(V::Lar.Points, EVs::Array{Lar.Cells})

model = V,EVs
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
