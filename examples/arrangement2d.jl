using LinearAlgebraicRepresentation
using Plasm,SparseArrays
Lar = LinearAlgebraicRepresentation


# input generation
V,EV = Lar.randomcuboids(30, .35)
V = Plasm.normalize(V,flag=true)
Plasm.view(Plasm.numbering(.05)((V,[[[k] for k=1:size(V,2)], EV])))

# subdivision of input edges
W = convert(Lar.Points, V')
cop_EV = Lar.coboundary_0(EV::Lar.Cells)
cop_EW = convert(Lar.ChainOp, cop_EV)
V, copEV = Lar.planar_arrangement_1(W::Lar.Points, cop_EW::Lar.ChainOp)

# compute containment graph of components
bicon_comps = Lar.Arrangement.biconnected_components(copEV)
@show bicon_comps
@show V
@show Matrix(copEV)

# visualization of component graphs
EW = Lar.cop2lar(copEV)
W = convert(Lar.Points, V')
hpcs = [ Plasm.lar2hpc(W,[EW[e] for e in comp]) for comp in bicon_comps ]
Plasm.view([ Plasm.color(Plasm.colorkey[(k%12)==0 ? 12 : k%12])(hpcs[k]) for k=1:(length(hpcs)) ])

# computation of 2-cells and 2-boundary
V, copEV, copFE = Lar.planar_arrangement_2(V, copEV, bicon_comps)




# final visualization
triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
V = convert(Lar.Points, V')
FVs = convert(Array{Lar.Cells}, triangulated_faces)
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

EVs = Lar.FV2EVs(copEV, copFE) # polygonal face boundaries
EVs = convert(Array{Array{Array{Int64,1},1},1}, EVs)
Plasm.viewcolor(V::Lar.Points, EVs::Array{Lar.Cells})

model = V,EVs
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))









"""
	numbering(V, copEV::Lar.ChainOp, copFE::Lar.ChainOp)::Plasm.Hpc

Produce the numbered `Hpc` of `planar_arrangement()` 2D output. 
Vertices in `V` are stored by row.

# Example

```julia
using LinearAlgebraicRepresentation
using Plasm, SparseArrays
Lar = LinearAlgebraicRepresentation

V,EV = Lar.randomcuboids(7, 1.0);
V = Plasm.normalize(V,flag=true);
W = convert(Lar.Points, V');
cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
V, copEV, copFE = Lar.planar_arrangement(W, cop_EW);

Plasm.view( numbering(0.05)((V, copEV, copFE)) )
```
"""
function numbering(scaling=0.1)
	function numbering0((V, copEV, copFE))
		VV = [[k] for k=1:size(V,1)]
		EV = [findnz(copEV[h,:])[1] for h=1:size(copEV,1)]
		FV = [collect(Set(cat(EV[e] for e in findnz(copFE[i,:])[1]))) for i=1:size(copFE,1)]
		FV = convert(Array{Array{Int64,1},1}, FV)
		model = (convert(Lar.Points, V'), Lar.Cells[VV,EV,FV])
		return Plasm.numbering(scaling)(model)
	end
	return numbering0
end

