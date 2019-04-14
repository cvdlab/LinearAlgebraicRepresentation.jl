using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

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
    directions = mapslices(normalize, rand(n, 2) - .5*ones(n, 2), dims=2)
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

function rubiks_example(ncubes = 3)
    V = Float64[
        0 0 0; 0 1 0;
        1 1 0; 1 0 0;
        0 0 1; 0 1 1;
        1 1 1; 1 0 1
    ]

    EV = sparse(Int8[
        -1  1  0  0  0  0  0  0;
        0 -1  1  0  0  0  0  0;
        0  0 -1  1  0  0  0  0;
        -1  0  0  1  0  0  0  0;
        -1  0  0  0  1  0  0  0;
        0 -1  0  0  0  1  0  0;
        0  0 -1  0  0  0  1  0;
        0  0  0 -1  0  0  0  1;
        0  0  0  0 -1  1  0  0;
        0  0  0  0  0 -1  1  0;
        0  0  0  0  0  0 -1  1;
        0  0  0  0 -1  0  0  1;
    ])

    FE = sparse(Int8[
        1  1  1 -1  0  0  0  0  0  0  0  0;
        0  0  0  0  0  0  0  0 -1 -1 -1  1;
        -1  0  0  0  1 -1  0  0  1  0  0  0;
        0 -1  0  0  0  1 -1  0  0  1  0  0;
        0  0 -1  0  0  0  1 -1  0  0  1  0;
        0  0  0  1 -1  0  0  1  0  0  0 -1;
    ])

    cube = [Lar.Points(V), EV, FE]
    cubesRow = [Lar.Points(zeros(0,3)),spzeros(Int8,0,0),spzeros(Int8,0,0)]

    println(map(typeof, cube))
    println(map(typeof, cubesRow))


    for i in 1:ncubes
        cubesRow = Lar.skel_merge(cubesRow..., cube...)
        cube[1] = Lar.Points(cube[1] + [zeros(8) zeros(8) ones(8)])
    end

    cubesRow = collect(cubesRow)
    cubesPlane = cubesRow
    num = size(cubesRow[1], 1)
    for i in 1:ncubes
        cubesPlane = Lar.skel_merge(cubesPlane..., cubesRow...)
        cubesRow[1] = Lar.Points(cubesRow[1] + [zeros(num) ones(num) zeros(num)])
    end

    cubesPlane = collect(cubesPlane)
    cubesCube = cubesPlane
    num = size(cubesPlane[1], 1)
    for i in 1:ncubes
        cubesCube = Lar.skel_merge(cubesCube..., cubesPlane...)
        cubesPlane[1] = Lar.Points(cubesPlane[1] + [ones(num) zeros(num) zeros(num)])
    end

    println("Arranging a cube of ", ncubes^3," cubes...")
    rubik = Lar.spatial_arrangement(cubesCube...)
    println("DONE")

    rubik = Lar.Points(rubik[1] - (.5*ncubes)), rubik[2:3]...
    c = cos(pi/6); s = sin(pi/6)
    M1 = [1  0 0; 0 c -s; 0 s c]
    M2 = [c -s 0; s c  0; 0 0 1]
    rot_rubik = Lar.Points(rubik[1]*M1*M2), rubik[2:3]...

    println("Arranging two rubik cubes...")
    two_rubiks = Lar.skel_merge(rubik..., rot_rubik...)
    println("DONE")

    arranged_rubiks = Lar.spatial_arrangement(two_rubiks...)
end

