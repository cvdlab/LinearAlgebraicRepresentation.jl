
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
    edge_boundary = Vector{Bool}(undef,10)
    counters = zeros(size(edges_list,1))
    for e in triaedges
        counters[e]+=1
    end
    for e in triaedges
        if counters[e] == 1 # outer edge
            edge_boundary[e] = false
        elseif counters[e] == 2 # inner edge
            edge_boundary[e] = true
        else
            error("in computing inner loops")
        end
    end
    # compute inner triangles
    inneredges = Array{Array{Int64,1},1}()
    for (k,value) in enumerate(edge_boundary)
       if value==true
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

    # all triangles are generated counterclockwise by Triangle.jl
    # Hence, close points to the left of directed hole edges must
    # be inner to the hole, and hence external to polygon, so that
    # the corresponding triangle is cancelled from the triangle array.
    # This must be repeated for each triangle classified as a hole.

    holes_list = Array{Array{Float64,1},1}()
    epsilon = 1e-5
    for (u,v,w) in holes  # each hole is a triangle
        a = normalize(V[:,v] - V[:,u]) #    unit side vector of first edge
        b = [-a[2,:]; a[1,:]]   #   unit orthogonal vector
        m = (V[:,u] + V[:,v])/2 #   middle side point
        pointinhole = m + epsilon * m  # very close to hole's first edge
        push!(holes_list, pointinhole)
    end
    triangles = [tria for tria in triangles if !(tria in holes)]



points = convert(Array{Float64,2}, V')
points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
edges_list = convert(Array{Int64,2}, hcat(EV...)')
edge_boundary = [true for k=1:size(edges_list,1)]
holes_list = convert(Array{Float64,2}, hcat(holes_list...)')
#triangles = Triangle.constrained_triangulation(points,points_map,edges_list,edge_boundary,holes_list)
triangles = [tria for tria in triangles if !(tria in holes)]


Plasm.view(Plasm.numbering(0.35)((V,[[[k] for k=1:size(V,2)], ev_nodups, triangles])))
Plasm.view(V,triangles)






#
# points = [0. 0.; 4. 0.; 2. 3.; 8. 0.; 6. 3.; 4. 6.]
# points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
# edges_list = Array{Int64,2}([1 2; 2 3; 3 1; 2 4; 4 5; 5 2; 3 5; 5 6; 6 3])
# edge_boundary = [false,true,false,false,false,true,true,false,false]
# holes_list = [4. 2.]
# triangles = Triangle.constrained_triangulation(points,points_map,edges_list,edge_boundary,holes_list)
