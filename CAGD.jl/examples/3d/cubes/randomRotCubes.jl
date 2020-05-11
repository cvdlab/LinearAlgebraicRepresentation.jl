todisplay = VERSION <= VersionNumber("1.2") ? true : false

using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
using Random
if todisplay
    using ViewerGL
    GL = ViewerGL
	include("../views.jl")
end

ncubes = 3
maxVisual = 10

function generateRotCubes(n = 2; axis = 3, seed = 1234)

    function rot(axis)
        idx = [true, true, true]
        if axis == 1
            idx = .!idx
            idx[rand(1:3)] = true
        elseif axis == 2
            idx[rand(1:3)] = false
        elseif axis != 3
            throw(ArgumentError("rot must be 1, 2 or three"));
        end
        x = Lar.r(rand()*idx[1], 0, 0)
        y = Lar.r(0, rand()*idx[2], 0)
        z = Lar.r(0, 0, rand()*idx[3])
        return x*y*z
    end

    Random.seed!(seed)
    V,(_,EV,FV,_) = Lar.cuboid([1,1,1],true,[-1,-1,-1])
    cube = V,FV,EV

    carry = Lar.Struct( cat([[rot(axis), cube] for k=1:n]) )

    V, FV, EV = Lar.struct2lar(carry)

    if todisplay  GL.VIEW([ GL.GLGrid(V,FV,GL.COLORS[1],0.5) ])  end

    m = CAGD.Model(V)
    CAGD.addModelCells!(m, 1, EV, signed = true)
    cFE = convert(Lar.ChainOp, Lar.coboundary_1(m.G, FV, EV))
    CAGD.addModelCells!(m, 2, cFE)
    return m
end


model = generateRotCubes(ncubes)

#=
atol = 1e-9;
split_model = CAGD.pairwise_decomposition(model, atol = atol)
congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)
if todisplay  displayModel(congr_model)  end
=#

A_model = CAGD.spatial_arrangement(model, atol = 1e-6, outerCell = true)

if todisplay
    V,CVs,FVs,EVs = Lar.pols2tria(A_model.G, A_model.T[1], A_model.T[2], A_model.T[3]);
    GL.VIEW(GL.GLExplode(V, FVs,         1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  1.0,1.0,1.0,  99,1.0));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  1.0,1.0,1.0,  99,0.5));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  5.0,5.0,5.0,  99,0.5));
end

if todisplay
    for i = 2 : min(length(CVs), maxVisual)
        GL.VIEW([
            GL.GLExplode(V,CVs[i:i], 1.0,1.0,1.0, 2, 0.9)
            GL.GLExplode(V,CVs[1:1], 1.0,1.0,1.0, 1, 0.1)
            GL.GLFrame
        ])
    end
end

