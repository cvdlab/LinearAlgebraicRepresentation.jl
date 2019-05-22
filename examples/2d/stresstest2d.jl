using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
using IntervalTrees
using SparseArrays
using NearestNeighbors

function input_collection(data::Array)::Lar.LAR
	assembly = Lar.Struct(data)
	return Lar.struct2lar(assembly)
end

function stresstest2d()
	V,(_,EV,FV) = Lar.cuboidGrid([4,1],true);
	W,(_,EW,FW) = Lar.cuboidGrid([1,5],true);
	mycircle(r,n) = Lar.circle(r)(n)
	data2d1 = (V,EV)
	data2d2 = Lar.Struct([ Lar.t(2,2), Lar.r(pi/3), Lar.t(-1.5,-2.5), (W,EW) ])
	data2d3 = Lar.Struct([ Lar.t(2,2), mycircle(2.5,160) ])
	data2d4 = Lar.Struct([ Lar.t(3.5,3.5), mycircle(.25,160) ])
	data2d5a = Lar.Struct([ Lar.t(5,3.5), mycircle(.75,160) ])
	data2d5 = Lar.Struct([ Lar.t(5,3.5), mycircle(.5,160) ])
	data2d6 = Lar.Struct([ Lar.t(5,3.5), mycircle(.25,160) ])
	data = [ data2d1, data2d2, data2d3, data2d4, data2d5,  data2d5a, data2d6 ]
	model2d = input_collection( [ Lar.Struct(data), Lar.t(-pi/6,pi/3), Lar.Struct(data) ] )
	V,EV = model2d
	Plasm.view( V, EV )

	W = convert(Lar.Points, V')
	cop_EV = Lar.coboundary_0(EV::Lar.Cells)
	cop_EW = convert(Lar.ChainOp, cop_EV)
	V, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)

	triangulated_faces = Lar.triangulate2D(V, [copEV, copFE])
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	V = convert(Lar.Points, V')
	Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

	W, copEV, copFE = Lar.Arrangement.planar_arrangement(W::Lar.Points, cop_EW::Lar.ChainOp)
	EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
	V = convert(Lar.Points, W')
	Plasm.viewcolor(V::Lar.Points, EVs::Array{Lar.Cells})

	model = V,EVs
	Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
end

stresstest2d()
