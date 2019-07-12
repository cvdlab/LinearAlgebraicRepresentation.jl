using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL

V,(VV,EV,FV) = Lar.simplex(2, true)
triangle = (V,EV,FV)
model = Lar.Struct([ triangle, Lar.t(.15,.15), Lar.s(.5,.5), triangle ])
V,EV,FV = Lar.struct2lar(model)
VV = [[k] for k=1:size(V,2)]
GL.VIEW(GL.numbering(.4)((V,[VV, EV, FV]),GL.COLORS[1],0.2));

triangles = Lar.triangulate2d(V, EV)
GL.VIEW([ GL.GLGrid(V,triangles) ])

# generate edges
ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
# remove duplicated edges from triangulation
nodupev = collect(Set(ev))
GL.VIEW(GL.numbering(.4)((V,[VV, nodupev])));
