using LinearAlgebraicRepresentation
using SparseArrays
Lar = LinearAlgebraicRepresentation
L = Lar
using ViewerGL
GL = ViewerGL

function pols2tria(W, copEV, copFE, copCF) # W by columns
	V = convert(Lar.Points,W')
	triangulated_faces = Lar.triangulate(V, [copEV, copFE])
	EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
	FVs = convert(Array{Lar.Cells}, triangulated_faces)
	CVs = []
	for cell in 1:copCF.m
		obj = []
        for f in copCF[cell, :].nzind
            triangles = triangulated_faces[f]
			append!(obj, triangles)
        end
		push!(CVs,obj)
    end
	V = convert(Lar.Points,V')
	return V,CVs,FVs,EVs
end

W, copEV, copFE, copCF = twocubegrids(2,2,2);
V,CVs,FVs,EVs = pols2tria(W, copEV, copFE, copCF)

meshes = GL.GLExplode(V,CVs[2:end],3,3,3,1,0.2)
GL.VIEW(meshes);
