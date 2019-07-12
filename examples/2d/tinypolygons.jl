using LinearAlgebraicRepresentation, ViewerGL
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using QHull, LinearAlgebra, DataStructures


function makecycle(points,edges_list)
    edgedict = SortedDict(zip(edges_list[:,1],edges_list[:,2]))
    polygon = [1]
    while true
        h = polygon[end]
        push!(polygon,edgedict[h])
        if polygon[end] == polygon[1] break end
    end
    return polygon # as ordered sequence of vertex indices
end

function vectors(V,polygon)
    vect = Array{Array{Float64,1},1}(undef,length(polygon))
    vect[1] = [0.0,0.0]
    for k=2:length(polygon)
        vect[k] = (V[polygon[k],:] - V[polygon[k-1],:])
    end
    return vect
end

function paramsum(vect)
    out = Array{Float64}(undef,length(vect))
    out[1] = 0.0
    for k=2:length(vect)
        out[k] = out[k-1] + norm(vect[k])
    end
    return out ./ out[end]
end

# problematic dataset
points = [-0.709583 -0.632635; 0.520341 0.561664; 0.918992 1.04452; 0.571694 0.592953; -0.120997 -0.0175708; -0.595547 -0.530965; 0.277832 0.413903; -0.980521 -0.823759]
points_map = [1, 2, 3, 4, 5, 6, 7, 8]
edges_list = [1 6; 3 8; 7 2; 2 4; 4 3; 5 7; 6 5; 8 1]

# problematic polygon
npoints = size(points,1)
polygon = makecycle(points,edges_list)
vects = vectors(points,polygon)
params = paramsum(vects) 
vdict = Dict(zip([points[k,:] for k=1:npoints], 1:npoints))

# convex hull polygon
ch = QHull.chull(points)
chverts = ch.vertices
outverts = setdiff(1:npoints, chverts)
# place the outverts on the convex hupp polygon
