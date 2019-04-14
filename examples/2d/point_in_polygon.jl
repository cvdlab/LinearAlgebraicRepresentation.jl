using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm
using NearestNeighbors
using PyCall

p = PyCall.pyimport("pyplasm")

filename = "/Users/paoluzzi/Documents/dev/Plasm.jl/test/svg/tile.svg"
V,EV = Plasm.svg2lar(filename)
Plasm.view(Plasm.numbering(.2)((V,[[[k] for k=1:size(V,2)], EV])))

function pointsinout(V,EV, n=10000)
	result = [Plasm.lar2hpc(V,EV)]
	classify = Lar.pointInPolygonClassification(V,EV)
	Threads.@threads for k=1:n
		queryPoint = [rand(),rand()]    
		inOut = classify(queryPoint)
		# println("k = $k, queryPoint = $queryPoint, inOut = $inOut")
		if inOut=="p_in"
			push!(result, p["MK"](queryPoint))
		elseif inOut=="p_out"
			push!(result, Plasm.color("red")(p["MK"](queryPoint)))
		elseif inOut=="p_on"
			push!(result, Plasm.color("green")(p["MK"](queryPoint)))
		end
	end
	return result
end

result = @time pointsinout(V,EV);
Plasm.view(p["STRUCT"](result))


