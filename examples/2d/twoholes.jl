using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Triangle
V,(VV,EV,FV) = Lar.simplex(2, true)
W,(WW,EW,FW) = Lar.cuboid([1,1], true)
triangle = (V,EV)
square = (W,EW)
model = Lar.Struct([
            Lar.Struct([
                square,
                Lar.t(.175,.175), Lar.s(.5,.5),
                triangle
            ]),
            Lar.t(.825,.825), Lar.s(-.5,-.5),
            triangle
        ])
V,EV = Lar.struct2lar(model)
Plasm.view(Plasm.numbering(0.5)((V,[[[k] for k=1:size(V,2)], EV])))

points = convert(Array{Float64,2}, V')
points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
edges_list = convert(Array{Int64,2}, hcat(EV...)')
edge_boundary = [true for k=1:size(edges_list,1)]
#triangles = Triangle.constrained_triangulation(points,points_map,edges_list,edge_boundary)
triangles = Triangle.constrained_triangulation(points,points_map,edges_list)
ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))

Plasm.view(Plasm.numbering(0.5)((V,[[[k] for k=1:size(V,2)], ev_nodups])))

ev_nodups = collect(Set(ev))
edge_dict = Dict(zip(EV,1:length(EV)))
outeredges = [edge_dict[[u,v]] for (u,v) in ev if haskey(edge_dict, [u,v]) ]
originaledges = [edge_dict[[u,v]] for (u,v) in ev_nodups if haskey(edge_dict, [u,v]) ]

edge_boundary = Vector{Bool}(undef,10)
counters = zeros(size(edges_list,1))
for e in outeredges counters[e]+=1 end
for e in outeredges
    if counters[e] == 1
        edge_boundary[e] = true
    elseif counters[e] == 2
        edge_boundary[e] = false
    else
        error("in computing inner loops")
    end
end

edge_boundary

Plasm.view(Plasm.numbering(0.5)((V,[[[k] for k=1:size(V,2)], EV])))








points = [0. 0.; 4. 0.; 2. 3.; 8. 0.; 6. 3.; 4. 6.]
points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
edges_list = Array{Int64,2}([1 2; 2 3; 3 1; 2 4; 4 5; 5 2; 3 5; 5 6; 6 3])
edge_boundary = [false,true,false,false,false,true,true,false,false]
holes_list = [4. 2.]
triangles = Triangle.constrained_triangulation(points,points_map,edges_list,edge_boundary,holes_list)
