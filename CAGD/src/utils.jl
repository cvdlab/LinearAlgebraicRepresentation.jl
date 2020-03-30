function build_projection_matrix(vs::Lar.Points; y0 = false)::Matrix{Float64}
    u1 = vs[:, 2] - vs[:, 1]
    i = 3
    u2 = vs[:, i] - vs[:, 1]
    u3 = Lar.cross(u1, u2)
    while u3 == [0.0, 0.0, 0.0]
        i += 1
        u2 = vs[:, i] - vs[:, 1]
        u3 = Lar.cross(u1, u2)
    end
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
    # Faces are sorted by the `angles` their xy-ray form w.r.t. the z = 0 positive ray
    angles = zeros(length(faces))
    for fidx = 1 : length(faces)
        # Retrieve ray orientation either on xy or xz.
        # xy is chosen unless the face is normal to y axis.
        oraxis = sum(1 .- isapprox.(PG[2, FVs[fidx]], 0.0)) > 0 ? 2 : 3
        # Search for a consisten `v`
        for v in FVs[fidx]
            if !isapprox(PG[2:3, v], [0.,0.])
                # Check whether face is on the side described by `v`'s ray or not
                # `side_pt` is located near the middle of `τ` on the side of `v`
                side_pt = [0.5 sign(PG[oraxis, v])*1e-10]
                # Take the edges of the face as array of arrays
                EVs = EV[model.T[2][faces[fidx], :].nzind]
                other_side=Lar.pointInPolygonClassification(PG[[1,oraxis],:], EVs)(side_pt)=="p_out"
                # Copute the angle of `v`'s ray and reorient it on the face side
                angles[fidx] = (-1)^other_side * atan(PG[3, v], PG[2, v])
                break
            end
        end
    end

    ord_faces = sortperm(angles)
    
    return faces[ord_faces]

    #=
    # Build angles for each petal that is not spanned by [0. 0. 1.] (i.e. τ)
    # and relates a face to a representative angle (all should be equal)
    angles = [
        [
            faces[idx],
            [
                (-1)^(0)*
        #   mod(________________________, 2π)  if [0, 2π] instead of [-π, +π]
                atan(PG[3, v], PG[2, v])
                for v in FVs[idx]
                if !isapprox(PG[2:3, v], [0.,0.])
            ][1]
        ]
        for idx = 1 : length(faces)
    ]
    # Sort the angles (right hand order)
    sort!(angles, lt = (a, b) -> a[2] > b[2])
    # return the petals idxs only
    return map(a->Int(a[1]), angles)
    =#
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