using LinearAlgebraicRepresentation
using Plasm


V = [[406.177, -437.882], [341.471, -352.0], [208.529, -352.0], [114.412, -420.945], [297.64, -433.177], [317.14, -456.319], [273.14, -447.445], [304.999, -406.118], [382.647, -336.706], [288.529, -466.117], [335.588, -500.234], [370.883, -297.882], [349.705, -324.941], [362.059, -406.118], [317.14, -336.706], [190.882, -297.882], [341.471, -290.824], [183.823, -406.118], [333.442, -389.239], [161.471, -460.236], [250.882, -397.882], [239.706, -352.0], [349.705, -426.782], [297.64, -330.824], [321.471, -317.294], [288.529, -247.294], [297.64, -386.473], [226.176, -452.0], [240.294, -534.353], [341.471, -473.445]]

EV = [(18, 22), (4, 18), (0, 10), (23, 26), (12, 16), (3, 17), (16, 24), (6, 7), (6, 29), (20, 21), (19, 20), (9, 28), (13, 29), (5, 22), (27, 28), (7, 14), (17, 19), (15, 25), (4, 5), (8, 11), (21, 23), (9, 10), (2, 3), (11, 12), (0, 1), (26, 27), (13, 14), (24, 25), (1, 8), (2, 15)]


W = hcat(V...)
EW = [[v1,v2]+1 for (v1,v2) in EV]
w1,w2 = minimum(W[1,:]), minimum(W[2,:])
W = hcat([W[:,k]-[w1;w2]  for k=1:size(W,2)]...)
sx = maximum(W[1,:])
sy = maximum(W[2,:])
V,EV = Plasm.apply(Lar.s(1/sx,1/sy))((W,EW))
Plasm.view(V,EV)


using PyCall
@pyimport pyplasm as p

result = PyObject[]
classify = pointInPolygonClassification(V,EV)
p_in = p_out = p_on = 0
for k=1:10000
    queryPoint = [rand(),rand()]
    inOut = classify(queryPoint)
    if inOut=="p_in" 
    	p_in += 1
    	push!(result, p.COLOR(p.WHITE)(p.MK(queryPoint)) )
    elseif inOut=="p_out" 
    	p_out += 1
    	push!(result, p.COLOR(p.RED)(p.MK(queryPoint)) )
    elseif inOut=="p_on" 
    	p_on += 1
    	push!(result, p.COLOR(p.YELLOW)(p.MK(queryPoint)) )
    end
end
p.VIEW(p.STRUCT(result))







