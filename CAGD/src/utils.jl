function build_projection_matrix(vs::Lar.Points; y0 = false)::Matrix{Float64}

    u1 = vs[:, 2] - vs[:, 1]
    i = 2 + argmax(map(i -> norm(Lar.cross(u1, vs[:, i] - vs[:, 1])), 3 : size(vs, 2)))
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
    eval_ord_angle(model::CAGD.Model, dim::Int, lo_cell::Int)::Array{Int,1}

Evaluates the `dim`-cell circular ordering w.r.t the `dim-1` cell `lo_cell`.
"""
function eval_ord_angle(model, dim, lo_cell)
    if dim == 1
        return CAGD.eval_ord_edges(model, lo_cell)
    elseif dim == 2
        return CAGD.eval_ord_faces(model, lo_cell)
    end
    throw(ArgumentError("Dimension $dim not coded for ordering"))
end

function eval_ord_faces(model, τ)
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
    M = CAGD.build_projection_matrix(Gface, y0 = true)
    # Project the other points such that each face is normal to x axis
    PG = M[1:3, :] * [model.G; ones(1, size(model, 0, 2))]
    # Each face describe in yz a straight line.
    # Faces are sorted by the `angles` their xy-ray form w.r.t. the y = 0 positive ray
    angles = zeros(length(faces))
    for fidx = 1 : length(faces)
        # If the face lies on xz plan, then the angle of the projection on yz is either -π/2 or +π/2
        if sum(1 .- isapprox.(PG[2, FVs[fidx]], 0.0)) == 0
            # Look for a non trivial point (z coord non null)
            for v in FVs[fidx]  if !isapprox(PG[3, v], 0.0)
                # Check whether face is on the z-halfplane described by `v`'s ray or not
                # `side_pt` is located "near" the middle of `τ` on the halfplane of `v`
                side_pt = [PG[1, edge[2]]/2 sign(PG[3, v])*1e-10]
                # Take the edges of the face as array of arrays for pointInPolygonClassification
                EVs = EV[model.T[2][faces[fidx], :].nzind]
                # Check wheter v is on the other side or the same halfplane
                other_side=Lar.pointInPolygonClassification(PG[[1,2],:], EVs)(side_pt)=="p_out"
                # the 90° angle is either on the side of `v` or not
                angles[fidx] = (-1)^other_side * sign(PG[3, v]) * π / 2
                break
            end  end
        else
            # If not so, the face is a non vertical segment if projected on yz plan
            # Search for a consisten `v`, that is y coord not null
            for v in FVs[fidx]  if !isapprox(PG[2, v], 0.0)
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
    ord_faces = sortperm(-angles)
    
    return faces[ord_faces]
end

function eval_ord_edges(model, τ; G = model.G)
    edges = model.T[1][:, τ].nzind
    EVs = [setdiff((model.T[1][e, :]).nzind, τ)[1] for e in edges]
    angles = [
        [edges[idx] atan(G[2, EVs[idx]], G[1, EVs[idx]])]
        for idx = 1 : length(edges)
        if !isapprox(G[:, EVs[idx]], [0.,0.])
    ]
    sort!(angles, lt = (a, b) -> a[2] < b[2])
    return map(a->Int(a[1]), angles)
end