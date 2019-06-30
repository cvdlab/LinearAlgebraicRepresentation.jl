using JLD2
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
#using SparseArrays, LinearAlgebra
using ViewerGL
GL = ViewerGL

# created by surface_extraction_parallel_ircad01.jl
@load "test/liver01.jld2" V FV
FW = [[v1,v2,v4,v3] for (v1,v2,v3,v4) in FV]
GL.VIEW([ GL.GLGrid(V,FW,GL.COLORS[2]) ]);

# created by surface_extraction_parallel_ircad01.jl
@load "test/liver02.jld2" V FV
TV = cat([[[v1,v2,v3],[v1,v4,v3]] for (v1,v2,v3,v4) in FV])
TV = convert(Lar.Cells,TV)
GL.VIEW([ GL.GLGrid(V,TV,GL.COLORS[1]) ]);

Lar.lar2obj
