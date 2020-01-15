using LinearAlgebra
Lar = LinearAlgebraicRepresentation

function submanifold_mapping(vs)
    u1 = vs[2,:] - vs[1,:]
    u2 = vs[3,:] - vs[1,:]
    u3 = cross(u1, u2)
    T = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    T[4, 1:3] = - vs[1,:]
    M = Matrix{Float64}(LinearAlgebra.I, 4, 4)
    M[1:3, 1:3] = [u1 u2 u3]
    return T*M
end

function spatial_index(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)
    d = 3
    faces_num = size(FE, 1)
    IntervalsType = IntervalValue{Float64, Int64}
    boxes1D = Array{IntervalsType, 2}(undef, 0, d)

    for fi in 1:faces_num
        vidxs = (abs.(FE[fi:fi,:])*abs.(EV))[1,:].nzind
        intervals = map((l,u)->IntervalsType(l,u,fi),
        	Lar.bbox(V[vidxs, :])...)
        boxes1D = vcat(boxes1D, intervals)
    end
    trees = mapslices(IntervalTree{Float64, IntervalsType}, sort(boxes1D; dims=1), dims=1)

    function intersect_intervals(intervals)
        cells = Array{Int64,1}[]
        for axis in 1:d
            vs = map(i->i.value, intersect(trees[axis], intervals[axis]))
            push!(cells, vs)
        end
        mapreduce(x->x, intersect, cells)
    end

    mapping = Array{Int64,1}[]
    for fi in 1:faces_num
        cell_indexes = setdiff(intersect_intervals(boxes1D[fi, :]), [fi])
        push!(mapping, cell_indexes)
    end

    mapping
end

function face_int(V::Lar.Points, EV::Lar.ChainOp, face::Lar.Cell)
    vs = Lar.buildFV(EV, face)
    retV = Lar.Points(undef, 0, 3)

    visited_verts = []
    for i in 1:length(vs)
        o = V[vs[i],:]
        j = i < length(vs) ? i+1 : 1
        d = V[vs[j],:] - o

        # err = 10e-8
        err = 10e-4
        println("approximation error =$err on edge decomposition")
        if !(-err < d[3] < err)

            alpha = -o[3] / d[3]

            if -err <= alpha <= 1+err
                p = o + alpha*d

                if -err < alpha < err || 1-err < alpha < 1+err
                    if !(Lar.vin(p, visited_verts))
                        push!(visited_verts, p)
                        retV = [retV; reshape(p, 1, 3)]
                    end
                else
                    retV = [retV; reshape(p, 1, 3)]
                end
            end
        end

    end

    vnum = size(retV, 1)


    if vnum == 1
        vnum = 0
        retV = Lar.Points(undef, 0, 3)
    end
    enum = (÷)(vnum, 2)
    retEV = spzeros(Int8, enum, vnum)

    for i in 1:enum
        retEV[i, 2*i-1:2*i] = [-1, 1]
    end

    retV, retEV
end
