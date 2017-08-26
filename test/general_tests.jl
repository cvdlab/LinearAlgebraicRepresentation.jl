using Base.Test

function generate_perpendicular_lines(steps::Int, minlen, maxlen)
    V = zeros(0,2)

    function rec(o, d, s)
        if s == 0 return end

        a = (maxlen-minlen)*rand() + minlen
        p = o + a*d
        V = [V; o; p]

        b = (a-minlen)*rand() + minlen
        p = o + b*d
        rec(p, d, s-1)

        b = (a-minlen)*rand() + minlen
        p = o + b*d
        rec(p, perpendicular(d), s-1)
    end

    function perpendicular(vec)
        v = zeros(size(vec))
        v[1] = vec[2]
        v[2] = vec[1]
        return v
    end

    rec([0 0], [1 0], steps)
    rec([0 0], [0 1], steps)
    vnum = size(V, 1)
    enum = vnum >> 1
    EV = spzeros(Int8, enum, vnum)
    for i in 1:enum
        EV[i, i*2-1:i*2] = 1
    end
    V, EV
end


function generate_random_lines(n, points_range, alphas_range)
    origins = points_range[1] + (points_range[2]-points_range[1])*rand(n, 2)
    directions = mapslices(normalize, rand(n, 2) - .5*ones(n, 2), 2)
    alphas = alphas_range[1] + (alphas_range[2]-alphas_range[1])*rand(n)
    new_points = Array{Float64, 2}(n, 2)
    for i in 1:n
        new_points[i, :] = origins[i, :] + alphas[i]*directions[i, :]
    end
    V = [origins; new_points]
    EV = spzeros(Int8, n, n*2)
    for i in 1:n
        EV[i, i] = 1
        EV[i, n+i] = 1
    end
    V, EV
end

