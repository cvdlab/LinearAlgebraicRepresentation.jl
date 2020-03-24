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
    get_ord_angle(model::CAGD.Model, dim::Int, lo_cell::Int)::Array{Int,1}

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
#    faces = model.T[2][:, τ].nzind;
#    Gface, Gfaceidx = CAGD.getModelCellGeometry(model, 2, faces[1], true)
#    M = CAGD.build_projection_matrix(Gface, y0 = true)
#    PG = (M * [model.G; ones(1, size(model, 0, 2))])[1:2, :]
#    angles = [[σ, mod(atan(PG[σ, 2], PG[σ, 1]), 2π)] for σ in faces]
#    sort!(angles, lt = (a, b) -> a[2] < b[2])
#    return map(a->a[1], angles)

    # Get edges of τ's petals
    faces = model.T[2][:, τ].nzind
    # Get vertices idx of τ
    edge = model.T[1][τ, :].nzind
    # Get other vertices idx of petals (but τ's)
    FVs = [setdiff((model.T[1]'*model.T[2][f, :]).nzind, edge) for f in faces]
    # Get geometry of the first petal (with τ's vertices first)
    Gface = model.G[:, [edge..., FVs[1]...]]
    # Get transformation matrix that maps τ to (0,0,0)->(x,0,0)
    M = CAGD.build_projection_matrix(Gface, y0 = true)
    # Project the other points, taking (y, z) coordinates
    PG = (M * [model.G; ones(1, size(model, 0, 2))])[2:3, :]
    # Build angles for each petal
    angles = [[mod(atan(PG[2, v], PG[1, v]), 2π) for v in fv] for fv in FVs]
    # Purge angles equal to 0 (if a segment is collinear with τ)
    angles = [[α for α in a if !(α ≈ 0)] for a in angles]
    # Take a representative angle for each petal and build a map (petal -> angle)
    angles = [[faces[i], angles[i][1]] for i = 1 : length(faces)]
    # Sort the angles (right hand order)
    sort!(angles, lt = (a, b) -> a[2] > b[2])
    # return the petals idxs only
    return map(a->Int(a[1]), angles)
end