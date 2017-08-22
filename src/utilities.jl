const Verts = Array{Float64, 2}
const Cells = SparseMatrixCSC{Int8, Int}
const Cell = SparseVector{Int8, Int}

function bbox(vertices::Verts)
    minimum = mapslices(x->min(x...), vertices, 1)
    maximum = mapslices(x->max(x...), vertices, 1)
    minimum, maximum
end

function bbox_contains(container, contained)
    b1_min, b1_max = container
    b2_min, b2_max = contained
    all(map((i,j,k,l)->i<=j<=k<=l, b1_min, b2_min, b2_max, b1_max))
end
function face_area(V::Verts, EV::Cells, face::Cell)
    function triangle_area(triangle_points::Verts)
        ret = ones(3,3)
        ret[:, 1:2] = triangle_points
        return .5*det(ret)
    end

    area = 0
    ps = [0, 0, 0]

    for i in face.nzind
        edge = face[i]*EV[i, :]
        skip = false

        for e in edge.nzind
            if e != ps[1]
                if edge[e] < 0
                    if ps[1] == 0
                        ps[1] = e
                        skip = true
                    else
                        ps[2] = e
                    end
                else
                    ps[3] = e
                end
            else
                skip = true
                break
            end
        end

        if !skip
            area += triangle_area(V[ps, :])
        end
    end

    return area
end
function skel_merge(V1::Verts, V2::Verts, EV1::Cells, EV2::Cells)
    V = [V1; V2]
    EV = spzeros(Int8, EV1.m + EV2.m, EV1.n + EV2.n)
    EV[1:EV1.m, 1:EV1.n] = EV1
    EV[EV1.m+1:end, EV1.n+1:end] = EV2
    V, EV
end

