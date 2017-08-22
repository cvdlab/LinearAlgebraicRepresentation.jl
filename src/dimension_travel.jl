typealias Verts Array{Float64, 2}
typealias Cells SparseMatrixCSC{Int8, Int}
typealias Cell SparseVector{Int8, Int}
using IntervalTrees

function normalizer(V::Verts)
    d = size(V, 2)
    upper = mapslices(x->max(x...), V, 1)
    lower = mapslices(x->min(x...), V, 1)
    diff = upper-lower

    T = eye(d+1)
    T[d+1, 1:d] = -lower

    S = eye(d+1)
    S[1:d, 1:d] = diagm(vec(map(inv, diff)))

    mat = T*S
    mat, inv(mat)
end
function submanifold_mapping(V, face)
    p1, p2, p3 = map(i->V[face.nzind[i], :], 1:3)
    u1 = p2-p1
    u2 = p3-p1
    u3 = cross(u1, u2)
    T = eye(4)
    T[4, 1:3] = -p1
    M = eye(4)
    M[1:3, 1:3] = [u1 u2 u3]
    return T*M
end
function spatial_index(V::Verts, CV::Cells)
    d = size(V,2)
    cell_num = size(CV, 1)
    IntervalsType = IntervalValue{Float64, Int64}
    boxes1D = Array{IntervalsType, 2}(0, d)
    for ci in 1:cell_num
        intervals = map((l,u)->IntervalsType(l,u,ci), bbox(V,CV[ci, :])...)
        boxes1D = vcat(boxes1D, intervals)
    end
    trees = mapslices(IntervalTree{Float64, IntervalsType}, sort(boxes1D, 1), 1)
    
    function intersect_intervals(intervals)
        cells = Array{Int64,1}[]
        for axis in 1:d
            vs = map(i->i.value, intersect(trees[axis], intervals[axis]))
            push!(cells, vs)
        end
        mapreduce(x->x, intersect, cells)
    end
    
    mapping = Array{Int64,1}[]
    for ci in 1:cell_num
        cell_indexes = setdiff(intersect_intervals(boxes1D[ci, :]), [ci])
        push!(mapping, cell_indexes)
    end
    
    mapping
end

