"""
    bbox(vertices::Points)

The axis aligned bounding box of the provided set of n-dim `vertices`.

The box is returned as the couple of `Points` of the two opposite corners of the box.
"""
function bbox(vertices::Points)
    minimum = mapslices(x->min(x...), vertices, 1)
    maximum = mapslices(x->max(x...), vertices, 1)
    minimum, maximum
end

"""
    bbox_contains(container, contained)

Check if the axis aligned bounding box `container` contains `contained`.

Each input box must be passed as the couple of `Points` standing on the opposite corners of the box.
"""
function bbox_contains(container, contained)
    b1_min, b1_max = container
    b2_min, b2_max = contained
    all(map((i,j,k,l)->i<=j<=k<=l, b1_min, b2_min, b2_max, b1_max))
end

"""
    face_area(V::Points, EV::Cells, face::Cell)

The area of `face` given a geometry `V` and an edge topology `EV`.
"""
function face_area(V::Points, EV::Cells, face::Cell)
    return face_area(V, buildEV(EV), face)
end

function face_area(V::Points, EV::ChainOp, face::Cell)
    function triangle_area(triangle_points::Points)
        ret = ones(3,3)
        ret[:, 1:2] = triangle_points
        return .5*det(ret)
    end

    area = 0

    fv = buildFV(EV, face)

    verts_num = length(fv)
    v1 = fv[1]

    for i in 2:(verts_num-1)

        v2 = fv[i]
        v3 = fv[i+1]

        area += triangle_area(V[[v1, v2, v3], :])
    end

    return area
end

"""
    skel_merge(V1::Points, EV1::ChainOp, V2::VePointsrts, EV2::ChainOp)

Merge two **1-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, V2::Points, EV2::ChainOp)
    V = [V1; V2]
    EV = spzeros(Int8, EV1.m + EV2.m, EV1.n + EV2.n)
    EV[1:EV1.m, 1:EV1.n] = EV1
    EV[EV1.m+1:end, EV1.n+1:end] = EV2
    V, EV
end

"""
    skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)

Merge two **2-skeletons**
"""
function skel_merge(V1::Points, EV1::ChainOp, FE1::ChainOp, V2::Points, EV2::ChainOp, FE2::ChainOp)
    FE = spzeros(Int8, FE1.m + FE2.m, FE1.n + FE2.n)
    FE[1:FE1.m, 1:FE1.n] = FE1
    FE[FE1.m+1:end, FE1.n+1:end] = FE2
    V, EV = skel_merge(V1, EV1, V2, EV2)
    V, EV, FE
end

"""
    delete_edges(todel, V::Points, EV::ChainOp)

Delete edges and remove unused vertices from a **2-skeleton**.

Loop over the `todel` edge index list and remove the marked edges from `EV`.
The vertices in `V` which remained unconnected after the edge deletion are deleted too.
"""

function delete_edges(todel, V::Points, EV::ChainOp)
    tokeep = setdiff(collect(1:EV.m), todel)
    EV = EV[tokeep, :]
    
    vertinds = 1:EV.n
    todel = Array{Int64, 1}()
    for i in vertinds
        if length(EV[:, i].nzind) == 0
            push!(todel, i)
        end
    end

    tokeep = setdiff(vertinds, todel)
    EV = EV[:, tokeep]
    V = V[tokeep, :]

    return V, EV
end

function buildFV(EV::Cells, face::Cell)
    return buildFV(buildEV(EV), face)
end

function buildFV(EV::ChainOp, face::Cell)
    startv = -1
    nextv = 0
    edge = 0

    vs = Array{Int64, 1}()

    while startv != nextv
        if startv < 0
            edge = face.nzind[1]
            startv = EV[edge,:].nzind[face[edge] < 0 ? 2 : 1]
            push!(vs, startv)
        else
            edge = setdiff(intersect(face.nzind, EV[:, nextv].nzind), edge)[1]
        end
        nextv = EV[edge,:].nzind[face[edge] < 0 ? 1 : 2]
        push!(vs, nextv)

    end

    return vs[1:end-1]
end

function buildFV(EV::ChainOp, face)
    startv = face[1]
    nextv = startv

    vs = []
    visited_edges = []

    while true
        curv = nextv
        push!(vs, curv)

        edge = 0
        for edge in EV[:, curv].nzind
            nextv = setdiff(EV[edge, :].nzind, curv)[1]
            if nextv in face && (nextv == startv || !(nextv in vs)) && !(edge in visited_edges)
                break
            end
        end

        push!(visited_edges, edge)

        if nextv == startv
            break
        end
    end

    return vs
end

function buildFE(FV, edges)
    faces = []

    for face in FV
        f = []
        for (i,v) in enumerate(face)
            edge = [v, face[i==length(face)?1:i+1]]
            ord_edge = sort(edge)

            edge_idx = findfirst(e->e==ord_edge, edges)

            push!(f, (edge_idx, sign(edge[2]-edge[1])))
        end
        
        push!(faces, f)
    end

    FE = spzeros(Int8, length(faces), length(edges))

    for (i,f) in enumerate(faces)
        for e in f
            FE[i, e[1]] = e[2]
        end
    end

    return FE
end

function buildEV(edges, signed=true)
    setValue = [-1, 1]
    if signed == false
        setValue = [1, 1]
    end

    maxv = max(map(x->max(x...), edges)...)
    EV = spzeros(Int8, length(edges), maxv)

    for (i,e) in enumerate(edges)
        e = sort(collect(e))
        EV[i, e] = setValue
    end

    return EV
end





function build_bounds(edges, faces)
    EV = buildEV(edges)
    FV = map(x->buildFV(EV,x), faces)
    FE = buildFE(FV, edges)

    return EV, FE
end

function vin(vertex, vertices_set)
    for v in vertices_set
        if vequals(vertex, v)
            return true
        end
    end
    return false
end

function vequals(v1, v2)
    err = 10e-8
    return length(v1) == length(v2) && all(map((x1, x2)->-err < x1-x2 < err, v1, v2))
end

function triangulate(V::Points, EV::ChainOp, FE::ChainOp)

    triangulated_faces = Array{Any, 1}(FE.m)

    for f in 1:FE.m
        if f % 10 == 0
            print(".")
        end
        
        edges_idxs = FE[f, :].nzind
        edge_num = length(edges_idxs)
        edges = zeros(Int64, edge_num, 2)

        
        fv = buildFV(EV, FE[f, :])

        vs = V[fv, :]

        v1 = normalize(vs[2, :] - vs[1, :])
        v2 = [0 0 0]
        v3 = [0 0 0]
        err = 1e-8
        i = 3
        while -err < norm(v3) < err
            v2 = normalize(vs[i, :] - vs[1, :])
            v3 = cross(v1, v2)
            i = i + 1
        end
        M = reshape([v1; v2; v3], 3, 3)

        vs = (vs*M)[:, 1:2]
        
        for i in 1:length(fv)
            edges[i, 1] = fv[i]
            edges[i, 2] = i == length(fv) ? fv[1] : fv[i+1]
        end
        
        triangulated_faces[f] = TRIANGLE.constrained_triangulation(vs, fv, edges, fill(true, edge_num))

        tV = (V*M)[:, 1:2]
        
        area = face_area(tV, EV, FE[f, :])
        if area < 0 
            for i in 1:length(triangulated_faces[f])
                triangulated_faces[f][i] = triangulated_faces[f][i][end:-1:1]
            end
        end
    end

    return triangulated_faces
end


function point_in_face(origin, V::Points, ev::ChainOp)

    function pointInPolygonClassification(V,EV)

        function crossingTest(new, old, status, count)
        if status == 0
            status = new
            return status, (count + 0.5)
        else
            if status == old
                return 0, (count + 0.5)
            else
                return 0, (count - 0.5)
            end
        end
        end

        function setTile(box)
        tiles = [[9,1,5],[8,0,4],[10,2,6]]
        b1,b2,b3,b4 = box
        function tileCode(point)
            x,y = point
            code = 0
            if y>b1 code=code|1 end
            if y<b2 code=code|2 end
            if x>b3 code=code|4 end
            if x<b4 code=code|8 end
            return code
        end
        return tileCode
        end

        function pointInPolygonClassification0(pnt)
            x,y = pnt
            xmin,xmax,ymin,ymax = x,x,y,y
            tilecode = setTile([ymax,ymin,xmax,xmin])
            count,status = 0,0

            for k in 1:EV.m
                edge = EV[k,:]
                p1, p2 = V[edge.nzind[1], :], V[edge.nzind[2], :]
                (x1,y1),(x2,y2) = p1,p2
                c1,c2 = tilecode(p1),tilecode(p2)
                c_edge, c_un, c_int = xor(c1, c2), c1|c2, c1&c2
                
                if (c_edge == 0) & (c_un == 0) return "p_on" 
                elseif (c_edge == 12) & (c_un == c_edge) return "p_on"
                elseif c_edge == 3
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 count += 1 end
                elseif c_edge == 15
                    x_int = ((y-y2)*(x1-x2)/(y1-y2))+x2 
                    if x_int > x count += 1
                    elseif x_int == x return "p_on" end
                elseif (c_edge == 13) & ((c1==4) | (c2==4))
                        status, count = crossingTest(1,2,status,count)
                elseif (c_edge == 14) & ((c1==4) | (c2==4))
                        status, count = crossingTest(2,1,status,count)
                elseif c_edge == 7 count += 1
                elseif c_edge == 11 count = count
                elseif c_edge == 1
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 
                        status, count = crossingTest(1,2,status,count) 
                    end
                elseif c_edge == 2
                    if c_int == 0 return "p_on"
                    elseif c_int == 4 
                        status, count = crossingTest(2,1,status,count) 
                    end
                elseif (c_edge == 4) & (c_un == c_edge) return "p_on"
                elseif (c_edge == 8) & (c_un == c_edge) return "p_on"
                elseif c_edge == 5
                    if (c1==0) | (c2==0) return "p_on"
                    else 
                        status, count = crossingTest(1,2,status,count) 
                    end
                elseif c_edge == 6
                    if (c1==0) | (c2==0) return "p_on"
                    else 
                        status, count = crossingTest(2,1,status,count) 
                    end
                elseif (c_edge == 9) & ((c1==0) | (c2==0)) return "p_on"
                elseif (c_edge == 10) & ((c1==0) | (c2==0)) return "p_on"
                end
            end
            
            if (round(count)%2)==1 
                return "p_in"
            else 
                return "p_out"
            end
        end
        return pointInPolygonClassification0
    end
    
    return pointInPolygonClassification(V, ev)(origin) == "p_in"
end

################
### Tri Output
################

"""

# Example

```julia
	julia> cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
	[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
	[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
	
	julia> cube_2 = LARLIB.Struct([LARLIB.t(0,0,0.5), LARLIB.r(0,0,pi/3), cube_1])
	
	julia> V,FV,EV = LARLIB.struct2lar(LARLIB.Struct([ cube_1, cube_2 ]))
	
	julia> V,bases,coboundaries = LARLIB.chaincomplex(V,FV,EV)
	
	julia> (EV, FV, CV), (cscEV, cscFE, cscCF) = bases,coboundaries

	julia> FV # bases[2]
	18-element Array{Array{Int64,1},1}:
	 [1, 3, 4, 6]            
	 [2, 3, 5, 6]            
	 [7, 8, 9, 10]           
	 [1, 2, 3, 7, 8]         
	 [4, 6, 9, 10, 11, 12]   
	 [5, 6, 11, 12]          
	 [1, 4, 7, 9]            
	 [2, 5, 11, 13]          
	 [2, 8, 10, 11, 13]      
	 [2, 3, 14, 15, 16]      
	 [11, 12, 13, 17]        
	 [11, 12, 13, 18, 19, 20]
	 [2, 3, 13, 17]          
	 [2, 13, 14, 18]         
	 [15, 16, 19, 20]        
	 [3, 6, 12, 15, 19]      
	 [3, 6, 12, 17]          
	 [14, 16, 18, 20]        

	julia> CV # bases[3]
	3-element Array{Array{Int64,1},1}:
	 [2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20]
	 [2, 3, 5, 6, 11, 12, 13, 17]                    
	 [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]    
	 
	julia> cscEV # coboundaries[1]
	34×20 SparseMatrixCSC{Int8,Int64} with 68 stored entries: ...

	julia> cscFE # coboundaries[2]
	18×34 SparseMatrixCSC{Int8,Int64} with 80 stored entries: ...
	
	julia> cscCF # coboundaries[3]
	4×18 SparseMatrixCSC{Int8,Int64} with 36 stored entries: ...
	
	objs = LARLIB.lar2obj(V'::LARLIB.Points, cscEV::LARLIB.ChainOp, 
			cscFE::LARLIB.ChainOp, cscCF::LARLIB.ChainOp)
			
	open("./two_cubes.obj", "w") do f
    	write(f, objs)
	end


```
"""
function lar2obj(V::Points, EV::ChainOp, FE::ChainOp, CF::ChainOp)
    obj = ""
    for v in 1:size(V, 1)
        obj = string(obj, "v ", round(V[v, 1], 6), " ", round(V[v, 2], 6), " ", 
        	round(V[v, 3], 6), "\n")
    end

    print("Triangulating")
    triangulated_faces = triangulate(V, EV, FE)
    println("DONE")

    for c in 1:CF.m
    obj = string(obj, "\ng cell", c, "\n")
    for f in CF[c, :].nzind
        triangles = triangulated_faces[f]
        for tri in triangles
            t = CF[c, f] > 0 ? tri : tri[end:-1:1]
            obj = string(obj, "f ", t[1], " ", t[2], " ", t[3], "\n")
        end
    end
end

    return obj
end
function obj2lar(path)
    fd = open(path, "r")
    vs = Array{Float64, 2}(0, 3)
    edges = Array{Array{Int, 1}, 1}()
    faces = Array{Array{Int, 1}, 1}()

    while (line = readline(fd)) != ""
        elems = split(line)
        if length(elems) > 0
            if elems[1] == "v"

                x = parse(Float64, elems[2])
                y = parse(Float64, elems[3])
                z = parse(Float64, elems[4])
                vs = [vs; x y z]

            elseif elems[1] == "f"
                v1 = parse(Int, elems[2])
                v2 = parse(Int, elems[3])
                v3 = parse(Int, elems[4])

                e1 = sort([v1, v2])
                e2 = sort([v2, v3])
                e3 = sort([v1, v3])

                if !(e1 in edges)
                    push!(edges, e1)
                end
                if !(e2 in edges)
                    push!(edges, e2)
                end
                if !(e3 in edges)
                    push!(edges, e3)
                end

                push!(faces, sort([v1, v2, v3]))
            end
        end
    end

    close(fd)
    vs, build_bounds(edges, faces)...
end