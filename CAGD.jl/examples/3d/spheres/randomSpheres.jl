todisplay = VERSION <= VersionNumber("1.2") ? true : false

using Random
using SparseArrays
using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
if todisplay
	using ViewerGL
	GL = ViewerGL
	include("../examples/3d/views.jl")
end


function build_copFE(FV::Lar.Cells, EV::Lar.Cells)
	copFE = Lar.u_coboundary_1(FV, EV) # unsigned
	faceedges = [findnz(copFE[f,:])[1] for f=1:size(copFE,1)]

	f_edgepairs = Array{Array{Int64,1}}[]
	for f=1:size(copFE,1)
		edgepairs = Array{Int64,1}[]
		for v in FV[f]
			push!(edgepairs, [e for e in faceedges[f] if v in EV[e]])
		end
		push!(f_edgepairs, edgepairs)
	end
	copFE.=0
	for f=1:size(copFE,1)
		facepairs = sort(sort.(f_edgepairs[f]))
		# setting a reference
		copFE[f, facepairs[1][1]] = 1
		for (e1,e2) in facepairs
			@assert copFE[f, e1] != 0
			if copFE[f, e2] == 0
				copFE[f, e2] =
					-((EV[e1][1]==EV[e2][1])||(EV[e1][2]==EV[e2][2])) * copFE[f,e1]+
					((EV[e1][1]==EV[e2][2])||(EV[e1][2]==EV[e2][1])) * copFE[f,e1]
			end
		end
	end
	return copFE
end

Random.seed!(1234);

scaling = 1.2
no_sphere = 10
expl = 2.0
pt1, pt2 = 5, 10
#pt1, pt2 = 10, 10

carry = Array{Tuple{Array{Float64,2},Array{Array{Int64,1},1}}}(undef, no_sphere)
for i = 1 : no_sphere
	carry[i] = Lar.struct2lar(Lar.Struct([
		Lar.t(rand(0.5:1e-2:5, 3)...),
		Lar.sphere(rand(1:1e-2:3) * scaling)([pt1, pt2])
	]))
end
V, FV = Lar.struct2lar(Lar.Struct(carry))
FV = sort(sort.(FV))

model = CAGD.Model(V);
function getFaceEdges(fV::Array{Int,1})
    return [fV[1], fV[2]], [fV[1], fV[3]], [fV[2], fV[3]]
end
EV = unique([(map(f -> getFaceEdges(f), FV)...)...])
FE = [[3*i+1, 3*i+2, 3*i+3] for i = 0 : length(FV)-1]
CAGD.addModelCells!(model, 1, EV, signed = true)
CAGD.addModelCells!(model, 2, build_copFE(FV, EV))
model = CAGD.mergeModelVertices(model, signed_merge=true)

if todisplay  displayModel(model)  end

split_model = CAGD.pairwise_decomposition(model)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)

if todisplay  displayModel(congr_model)  end

gift_model = deepcopy(congr_model)
FC, bicon_comps = CAGD.tgw(congr_model, 3)
CAGD.addModelCells!(gift_model, 3, convert(Lar.ChainOp, FC'))

if todisplay  viewExplode(gift_model)  end