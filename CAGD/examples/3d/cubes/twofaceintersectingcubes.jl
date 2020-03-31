display = VERSION < VersionNumber("1.2") ? true : false

using LinearAlgebra
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
CAGD = CAGD
if display
    using ViewerGL
    GL = ViewerGL
end

function displayModel(model)
    V = model.G
    EV = Lar.cop2lar(model.T[1])
    FE = Lar.cop2lar(model.T[2])
    CF = Lar.cop2lar(model.T[3])
    FV = Lar.cop2lar(map(x -> Int8(x/2), abs.(model.T[2]) * abs.(model.T[1])))
    CV = []
    triangulated_faces = Lar.triangulate(convert(Lar.Points, V'), [model.T[1], model.T[2]])
    FVs = convert(Array{Lar.Cells}, triangulated_faces)
    EVs = Lar.FV2EVs(model.T[1], model.T[2])

    GL.VIEW([ GL.GLPol(V,EV, GL.COLORS[1]) ]);
    GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99));
    GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
end

V = [
    0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0
    0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0
    0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
]
EV = [
    [1, 2], [1, 3], [2, 4], [3, 4],
    [1, 5], [2, 6], [3, 7], [4, 8],
    [5, 6], [5, 7], [6, 8], [7, 8]
]
FE = [
    [1, 2, 3, 4], [9, 10, 11, 12],
    [1, 5, 6, 9], [2, 5, 7, 10], [3, 6, 8, 11], [4, 7, 8, 12]
]

m1 = CAGD.Model(V)
CAGD.addModelCells!(m1, 1, EV, signed = true)
CAGD.addModelCells!(m1, 2, FE, signed = true)

m2 = deepcopy(m1)
m2.G .*= 3
m2.G .-= [1.0; 1.0; -0.5]

model = CAGD.uniteModels(m1, m2)

if display  displayModel(model)  end

split_model = CAGD.pairwise_decomposition(model)

if display  displayModel(split_model)  end

congr_model = CAGD.mergeModelVertices(split_model, signed_merge=true)

if display  displayModel(congr_model)  end

gift_model, bicon_comps = CAGD.tgw(congr_model, 3)

if display  displayModel(gift_model)  end

# Each face is percurred two times
sum(gift_model.T[3], dims = 1)
