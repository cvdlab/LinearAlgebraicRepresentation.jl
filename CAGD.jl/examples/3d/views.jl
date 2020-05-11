todisplay = VERSION <= VersionNumber("1.2") ? true : false

using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using CAGD
if todisplay
    using ViewerGL
    GL = ViewerGL
end

function viewExplode(model; complete = false)
    V,CVs,FVs,EVs = Lar.pols2tria(model.G, model.T[1], model.T[2], model.T[3])

    GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
    GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
    GL.VIEW(GL.GLExplode(V,CVs[2 - complete:end],5,5,5,99,0.5));
end

function displayModel(model, exp = 1.5)
    V = model.G
    EV = Lar.cop2lar(model.T[1])
    FE = Lar.cop2lar(model.T[2])
    CF = Lar.cop2lar(model.T[3])
    FV = Lar.cop2lar(map(x -> Int8(x/2), abs.(model.T[2]) * abs.(model.T[1])))
    CV = []
    triangulated_faces = Lar.triangulate(convert(Lar.Points, V'), [model.T[1], model.T[2]])
    FVs = convert(Array{Lar.Cells}, triangulated_faces)
    EVs = Lar.FV2EVs(model.T[1], model.T[2])

    GL.VIEW([
        GL.GLAxis( GL.Point3d(0,0,0),GL.Point4d(1,1,1,1) )
        GL.GLPol(V,EV, GL.COLORS[1])
    ]);
    GL.VIEW(GL.GLExplode(V,FVs,exp,exp,exp,99));
    GL.VIEW(GL.GLExplode(V,EVs,exp,exp,exp,99,1));
end


function viewNum(model, numsize = 0.25)
	V = model.G
	EV = Lar.cop2lar(model.T[1])
	FV = Lar.cop2lar(map(x -> floor(Int8, x / 2), abs.(model.T[2]) * abs.(model.T[1])))
	VV = [[k] for k=1:size(V,2)]
	GL.VIEW([
        GL.GLAxis( GL.Point3d(0,0,0),GL.Point3d(1,1,1) )
        GL.numbering(numsize)((V,[VV, EV, FV]))
    ])
end



function viewNumEdge(model::CAGD.Model, τ::Int, numsize=0.25)
    V  = model.G
    EV = Lar.cop2lar(model.T[1])
    VV = [[k] for k=1:size(V,2)]
    faces = model.T[2][:, τ].nzind
    EVidx = unique(vcat([model.T[2][f, :].nzind for f in faces]...))
    EVs = Lar.cop2lar(model.T[1][EVidx, :])
    FVs = Lar.cop2lar(map(x -> floor(Int8, x / 2), abs.(model.T[2][faces, :]) * abs.(model.T[1])))
    VV = [[k] for k in unique(vcat(EVs...))]
    GL.VIEW([
        GL.GLAxis( GL.Point3d(0,0,0),GL.Point3d(1,1,1) )
        GL.GLPol(V, EV, GL.COLORS[1])
        GL.numbering(numsize)((V, [VV, EVs, FVs]))
    ])
end

function viewNumFace(model::CAGD.Model, τ::Int, numsize=0.25)
    V  = model.G
    EV = Lar.cop2lar(model.T[1])
    VV = [[k] for k=1:size(V,2)]
    EVidx = model.T[2][τ, :].nzind
    EVs = Lar.cop2lar(model.T[1][EVidx, :])
    FVs = Lar.cop2lar(map(x -> floor(Int8, x / 2), abs.(model.T[2][τ:τ, :]) * abs.(model.T[1])))
    VV = [[k] for k in unique(vcat(EVs...))]
    GL.VIEW([
        GL.GLFrame,
        GL.GLPol(V, EV, GL.COLORS[1]),
        GL.numbering(numsize)((V, [VV, EVs, FVs]))
    ])
end
