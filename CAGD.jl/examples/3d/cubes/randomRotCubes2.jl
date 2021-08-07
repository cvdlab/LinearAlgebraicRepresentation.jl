todisplay = VERSION <= VersionNumber("1.2") ? true : false

using LinearAlgebra
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
using Random
if todisplay
    using ViewerGL
    GL = ViewerGL
end

# Data 0
n = 3;
rngseed = 1234;
maxVisual = 100;
atol = 1e-6;

# Data 1
n = 7;
rngseed = 1111;
maxVisual = 10;
atol = 1e-4;


# Data 2
n = 8;
rngseed = 4321;
maxVisual = 10;
atol = 1e-4;

function generateRotCubes(n = 3; seed = 1234)
    rng = MersenneTwister(seed);
    V,(_,EV,FV,_) = Lar.cuboid([1,1,1],true,[-1,-1,-1])
    cube = V,FV,EV
    cubos = []
    for k = 1:n
        angles = rand(rng, 3) .* 2π
        R = Lar.r(angles[1],0,0) * Lar.r(0,angles[2],0) * Lar.r(0,0,angles[3])
        push!(cubos, Lar.apply(R, cube))
    end
    V,FV,EV = Lar.struct2lar(Lar.Struct(cubos))
    if todisplay  GL.VIEW([ GL.GLGrid(V,FV,GL.COLORS[1],0.5) ])  end

    model = CAGD.Model(V);
    CAGD.addModelCells!(model, 1, EV, signed = true)
    cFE = convert(Lar.ChainOp, Lar.coboundary_1(model.G, FV, EV))
    CAGD.addModelCells!(model, 2, cFE)

    return model
end

#== UNCOMMENT for step-by-step pipeline ==#  
#==

model = generateRotCubes(n, seed = rngseed)

atol = 1e-6;

if todisplay
    GL.VIEW([
        GL.GLPol(model.G, Lar.cop2lar(model.T[1]), GL.COLORS[1])
        GL.GLFrame
    ]);
end

split_model = CAGD.pairwise_decomposition(model, atol = atol)

if todisplay
    GL.VIEW([
        GL.GLPol(model.G, Lar.cop2lar(model.T[1]), GL.COLORS[1])
        GL.GLFrame
    ]);
end

congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)

if todisplay
    GL.VIEW([
        GL.GLPol(congr_model.G, Lar.cop2lar(congr_model.T[1]), GL.COLORS[1])
        GL.GLFrame
    ]);
end

A_model = deepcopy(congr_model);
FC, bicon_comps = CAGD.tgw(congr_model, 3)
CAGD.addModelCells!(A_model, 3, convert(Lar.ChainOp, FC'))

if todisplay
    GL.VIEW([
        GL.GLPol(A_model.G, Lar.cop2lar(A_model.T[1]), GL.COLORS[1])
        GL.GLFrame
    ]);
end
==#

model = generateRotCubes(n, seed = rngseed)
A_model = CAGD.spatial_arrangement(model, atol = atol, outerCell = true)

if todisplay
    V,CVs,FVs,EVs = Lar.pols2tria(A_model.G, A_model.T[1], A_model.T[2], A_model.T[3]);
    GL.VIEW(GL.GLExplode(V, EVs,         1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, FVs,         1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  1.0,1.0,1.0,  99,1.0));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  1.0,1.0,1.0,  99,0.5));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  5.0,5.0,5.0,  99,0.5));
    # Takes the 3-cell(s) with more 2-cells
    idxs = SparseArrays.sparse(length.(CVs) .== max(length.(CVs[2:end])...)).nzind;
    GL.VIEW([
        GL.GLExplode(V,CVs[idxs], 1.0,1.0,1.0, 2, 0.9)
        GL.GLExplode(V,CVs[1:1], 1.0,1.0,1.0, 1, 0.1)
    ])
    GL.VIEW([
        GL.GLFrame
        GL.GLExplode(V, CVs[idxs], 1.0,1.0,1.0, 2, 0.9)
    ])
end

# Visualise first maxVisual 3-cells within the complete structure
if todisplay
    for i = 2 : min(length(CVs), maxVisual)
        GL.VIEW([
            GL.GLFrame
            GL.GLExplode(V,CVs[i:i], 1.0,1.0,1.0, 2, 0.9)
            GL.GLExplode(V,CVs[1:1], 1.0,1.0,1.0, 1, 0.1)
        ])
    end
end