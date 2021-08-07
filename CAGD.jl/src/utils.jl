using LinearAlgebra

function build_projection_matrix(vs::Lar.Points; y0 = false)::Matrix{Float64}

#    u1 = LinearAlgebra.normalize(vs[:, 2] - vs[:, 1])
    u1 = vs[:, 2] - vs[:, 1]
    i = 2 + argmax(map(i -> norm(Lar.cross(u1, vs[:, i] - vs[:, 1])), 3 : size(vs, 2)))
#    u2 = LinearAlgebra.normalize(vs[:, i] - vs[:, 1])
#    u3 = LinearAlgebra.normalize(Lar.cross(u1, u2))
#    if y0  u2 = LinearAlgebra.normalize(Lar.cross(u1, u3))  end
    u2 = vs[:, i] - vs[:, 1]
    u3 = Lar.cross(u1, u2)
    if y0  u2 = Lar.cross(u1, u3)  end
    T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    T[1:3, 4] = - vs[:,1]
    M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    M[1:3, 1:3] = [u1 u2 u3]
    return M' * T
end


"""
    build_rototranslation_matrix(vs::Lar.Points)::Matrix{Float64}

Generates the rototranslation matrix that maps `vs[1,2]` -> ([0 0 0; x 0 0]')
"""
function build_rototranslation_matrix(vs::Lar.Points)::Matrix{Float64}

    Ivs = [vs; ones(1, size(vs, 2))]
    T   = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    Rx  = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    Ry  = Matrix{Float64}(LinearAlgebra.I, 4, 4)

    # Translate first point to Origin
    T[1:3, 4] = - vs[:,1]
    Tvs = T * Ivs

    # Rotate along x-axis to xz-plan the 1->2 edge (y = 0)
    #   |  1   0   0  |       sin(θ) = Py (= θ[1])
    #   |  0   c  -s  |   =>
    #   |  0   s   c  |       cos(θ) = Pz (= θ[2])
    if !isapprox(norm(Tvs[2:3, 2]), 0.0, atol = 1e-10)
        θ = LinearAlgebra.normalize(Tvs[2:3, 2])
        Rx[2, 2] =  θ[2]
        Rx[3, 3] =  θ[2]
        Rx[2, 3] = -θ[1]
        Rx[3, 2] =  θ[1]
    end
    RTvs = (Rx * T) * Ivs

    # Rotate along y-axis to x-axis the 1->2 edge (y = 0 && z = 0)
    #   |  c   0   s  |       sin(σ) = Pz (= σ[2])
    #   |  0   1   0  |   =>
    #   | -s   0   c  |       cos(σ) = Px (= σ[1])
    if !isapprox(norm(RTvs[[1, 3], 2]), 0.0, atol = 1e-10)
        σ = LinearAlgebra.normalize(RTvs[[1, 3], 2])
        Ry[1, 1] =  σ[1]
        Ry[3, 3] =  σ[1]
        Ry[1, 3] =  σ[2]
        Ry[3, 1] = -σ[2]
    end

    return Ry * Rx * T
end


"""
    eval_ord_angle(model::CAGD.Model, dim::Int, lo_cell::Int)::Array{Int,1}

Evaluates the `dim`-cell circular ordering w.r.t the `dim-1` cell `lo_cell`.
"""
function eval_ord_angle(model, dim, lo_cell; atol = 1e-7)
    if dim == 1
        return CAGD.eval_ord_edges(model, lo_cell, atol = atol)
    elseif dim == 2
        return CAGD.eval_ord_faces(model, lo_cell, atol = atol)
    end
    throw(ArgumentError("Dimension $dim not coded for ordering"))
end

function eval_ord_faces(model, τ; atol = 1e-7)
    proj_matrix = false
    # Get EV array of arrays representation for Lar.pointInPolygonClassification
    EV = Lar.cop2lar(model.T[1])
    # Get faces of τ's petals
    faces = model.T[2][:, τ].nzind
    # Get vertices idx of τ
    edge = model.T[1][τ, :].nzind
    # Get other vertices idx of petals (but τ's)
    FVs = [setdiff((abs.(model.T[1])'*abs.(model.T[2][f, :])).nzind, edge) for f in faces]
    # Get geometry of the first petal (with τ's vertices first)
    Gface = model.G[:, [edge..., FVs[1]...]]
    # Get transformation matrix that maps τ to (0,0,0)->(x,0,0)
    if proj_matrix
        M = CAGD.build_projection_matrix(Gface, y0 = true)
    else
        M = CAGD.build_rototranslation_matrix(Gface)
    end
    # Project the other points such that each face is normal to x axis
    PG = M[1:3, :] * [model.G; ones(1, size(model, 0, 2))]
    # Each face describe in yz a straight line.
    # Faces are sorted by the `angles` their xy-ray form w.r.t. the y = 0 positive ray
    angles = zeros(length(faces))
    for fidx = 1 : length(faces)
        # If the face lies on xz plan, then the angle of the projection on yz is either -π/2 or +π/2
        if sum(1 .- isapprox.(PG[2, FVs[fidx]], 0.0, atol = atol)) == 0
            # Look for a non trivial point (z coord non null)
            for v in FVs[fidx]  if !isapprox(PG[3, v], 0.0, atol = atol)
                # Check whether face is on the z-halfplane described by `v`'s ray or not
                # `side_pt` is located "near" the middle of `τ` on the halfplane of `v`
                side_pt = [PG[1, edge[2]]/2 sign(PG[3, v])*1e-10]
                # Take the edges of the face as array of arrays for pointInPolygonClassification
                EVs = EV[model.T[2][faces[fidx], :].nzind]
                # Check wheter v is on the other side or the same halfplane
                other_side=Lar.pointInPolygonClassification(PG[[1,3],:], EVs)(side_pt)=="p_out"
                # the 90° angle is either on the side of `v` or not
                angles[fidx] = (-1)^other_side * sign(PG[3, v]) * π / 2
                break
            end  end
        else
            # If not so, the face is a non vertical segment if projected on yz plan
            # Search for a consisten `v`, that is y coord not null
            for v in FVs[fidx]  if !isapprox(PG[2, v], 0.0, atol = atol)
                # Check whether the face is on the halfplane described by `v`'s ray or not
                #  to do so, consider the xy projection of the face and check wether
                #  `side_pt` (located "near" the middle of `τ` on the halfplane of `v`)
                #  is in the face projection 
                side_pt = [PG[1, edge[2]]/2 PG[2, edge[2]]+sign(PG[2, v])*1e-10]
                # Take the edges of the face as array of arrays for pointInPolygonClassification
                EVs = EV[model.T[2][faces[fidx], :].nzind]
                # Check wheter v is on the other side or the same halfplane
                other_side=Lar.pointInPolygonClassification(PG[[1,2],:], EVs)(side_pt)=="p_out"
                # Copute the angle of `v`'s ray, reoriented according the face halfplane
                angles[fidx] = atan((-1)^other_side * PG[3, v], (-1)^other_side * PG[2, v])
                break
            end  end
        end
    end

    # faces are sorted according to right-hand rule
    if proj_matrix angles = -angles end
    ord_faces = sortperm(angles)
    
    return faces[ord_faces]
end

function eval_ord_edges(model, τ; G = model.G, atol = 1e-7)
    edges = model.T[1][:, τ].nzind
    EVs = [setdiff((model.T[1][e, :]).nzind, τ)[1] for e in edges]
    angles = [
        [edges[idx] atan(G[2, EVs[idx]], G[1, EVs[idx]])]
        for idx = 1 : length(edges)
        if !isapprox(G[:, EVs[idx]], [0.,0.], atol = atol)
    ]
    sort!(angles, lt = (a, b) -> a[2] < b[2])
    return map(a->Int(a[1]), angles)
end

"""
    pointInSolid(
        model::CAGD.Model,
        point::Array{Float64,1},
        cycle::Int,
        signed = false
    )::Int8

Computes if `point` is internal (`1`), external (`-1`) or onto (`0`) `cycle` solid

Checks the incidences of the `x`-ray from `point` with `model.T[3][cycle,:]`
boundary faces to determine wether it is internal, external or onto it.

By default it consider the cycle to be the internal boundary of a cell.
If `signed` is set `true` then it also consider the orientation of the cycle
(_i.e_ if the cycle is the external boundary). The model needs being oriented.
"""
function pointInSolid(model, point, cycle, signed = false; atol = 1e-8)

    function doBBoxIntersect(bbox) #dimension free
        flag = point[1] <= bbox[1, 2]
        for d = 2 : length(point)
            flag &= bbox[d, 1] <= point[d] <= bbox[d, 2]
        end
        return flag
    end

    if signed
        throw(ArgumentError("signed not coded yet"))
    end
    # Vertices of `cycle`
    vs = CAGD.getModelCellVertices(model, 3, cycle)
    # Check if `point` is a vertex
    for vert in vs  if isapprox(vs, vert)  return 0  end  end
    # Reshaping vertices
    vs = hcat(vs...)
    # Bounding box of vs
    minimum = mapslices(x->min(x...), vs, dims=2)
    maximum = mapslices(x->max(x...), vs, dims=2)
    # x-line vertices (outside of bbox)
    vmin = point
    # adding a little delta to avoid coplanarity
    vmax = [maximum[1]+1; point[2] + 2*atol; point[3] + 2*atol]

    # initialization of counter
    counter = 0
    face_cycles = CAGD.getModelLoCell(model, 3, cycle)
    face_bboxes = CAGD.getModelBoundingBoxes(model, 2, face_cycles, atol = atol*1e-2)
    for face_idx = 1 : length(face_cycles)  if doBBoxIntersect(face_bboxes[face_idx])
        # To check incidence on each face, project it to z=0
        Gface, Gfaceidx = CAGD.getModelCellGeometry(model, 2, face_idx, true)
        # Build the projection matrix
        M = CAGD.build_projection_matrix(Gface)
        # Project the points
        PG = (M[1:3, :] * [model.G; ones(1, size(model, 0, 2))])
        Pvmin = M[1:3, :] * [vmin; 1]
        Pvmax = M[1:3, :] * [vmax; 1]

        if isapprox(Pvmin[3], 0.0, atol = atol) && isapprox(Pvmax[3], 0.0, atol = atol)
            # if point is collinear with face, then it is either on face or not.
            # If not, then the face has even contribute and can be avoided
            if CAGD.point_in_face(Pvmin[1:2], PG[1:2, :], model.T[1])
                return 0
            end
        else
            # If point is not collinear, retrieve the x-ray point on face plan
            Pvplan = Pvmin[1:2] + (Pvmax[1:2] - Pvmin[1:2]).*((Pvmin[3])/(Pvmax[3]-Pvmin[3]))
            # If it is in, then it adds one to the counter
            if CAGD.point_in_face(Pvplan, PG[1:2, :], model.T[1])
                # check it is really in and not on the border
                @assert Lar.point_in_face(Pvplan, PG[1:2, :], model.T[1])
                counter += 1
            end
        end
    end  end
    return counter % 2 == 1
end


"""
    getExternalBoundary(model::CAGD.Model, deg::Int, cycles::Lar.ChainOp, comp::Array{Int,1})::Int

Evaluates the external cycle of the component described by `deg`-`cycles` `comp`

Consider a component made of the `comp` cycles among `deg`-cells and
evaluates wich of those is the external one (in `comp` ordering).

#TODO improve the function checking if the longest cycle is the outer
"""
function getExternalBoundary(model, deg, cycles, comp)
    return findmax([length(cycles[:, comp[i]].nzind) for i = 1 : length(comp)])[2]
end