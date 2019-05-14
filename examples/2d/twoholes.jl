using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Triangle, Plasm, LinearAlgebra

"""
    triangulate2d(V, EV)

Contrained Delaunay Triangulation of LAR model (V,EV).
Discovery and removal of holes from triangulation, by comparing
original and generated edges.
"""
function triangulate2d(V, EV)
    # data for Constrained Delaunay Triangulation (CDT)
    points = convert(Array{Float64,2}, V')
    points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
    edges_list = convert(Array{Int64,2}, hcat(EV...)')
    edge_boundary = [true for k=1:size(edges_list,1)]
    triangles = Triangle.constrained_triangulation(points,points_map,edges_list)
    # edges of the triangulation
    ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
    # remove duplicated edges from triangulation
    ev_nodups = collect(Set(ev))
    ##Plasm.view(Plasm.numbering(0.35)((V,[[[k] for k=1:size(V,2)], ev_nodups])))
    # dictionary o original edges
    edge_dict = Dict(zip(EV,1:length(EV)))
    triaedges = [edge_dict[[u,v]] for (u,v) in ev if haskey(edge_dict, [u,v]) ]
    # subdivide original edges between inner and outer
    edge_boundary = Vector{Bool}(undef,length(triaedges))
    counters = zeros(size(edges_list,1))
    for e in triaedges
        counters[e]+=1
    end
    edge_boundary = Vector{Bool}(undef,length(triaedges))
    for e in triaedges
        edge_boundary[e] = counters[e] == 1 ? false : true
    end
    # compute inner triangles
    inneredges = Array{Array{Int64,1},1}()
    for (k,value) in enumerate(counters)
       if value==2
           push!(inneredges, EV[k], reverse(EV[k]))
       end
    end
    # compute hole(s): wheater all triangle edges are inneredges
    holes = Array{Array{Int64,1},1}()
    for (k,(u,v,w)) in enumerate(triangles)
        triangle = [[u,v],[v,w],[w,u]]
        if setdiff(triangle,inneredges)==[]
            push!(holes, triangles[k])
        end
    end
    triangles = [tria for tria in triangles if !(tria in holes)]
    return triangles
end


# input of primitive shapes
V,(VV,EV,FV) = Lar.simplex(2, true)
W,(WW,EW,FW) = Lar.cuboid([1,1], true)
triangle = (V,EV)
square = (W,EW)
# hierarchical assembly
model = Lar.Struct([
            Lar.Struct([
                square,
                Lar.t(.175,.175), Lar.s(.5,.5),
                triangle
            ]),
            Lar.t(.825,.825), Lar.s(-.5,-.5),
            triangle
        ])
# visualization of generated wire-frame model
V,EV = Lar.struct2lar(model)
Plasm.view(Plasm.numbering(0.35)((V,[VV, EV])))

triangles = triangulate2d(V, EV)
Plasm.view(V,triangles)

ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
# remove duplicated edges from triangulation
ev_nodups = collect(Set(ev))
Plasm.view(Plasm.numbering(0.35)((V,[VV, ev_nodups, triangles])))
