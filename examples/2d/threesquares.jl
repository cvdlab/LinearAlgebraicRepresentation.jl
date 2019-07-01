using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Triangle, Plasm, LinearAlgebra

# input of primitive shapes
W,(WW,EW,FW) = Lar.cuboid([1,1], true)
square = (W,EW)
# hierarchical assembly
model = Lar.Struct([
            Lar.Struct([
                square,
                Lar.t(.175,.175), Lar.s(.25,.25),
                square
            ]),
            Lar.t(.825,.825), Lar.s(-.25,-.25),
            square
        ])
# visualization of generated wire-frame model
V,EV = Lar.struct2lar(model)
VV = [[k] for k=1:size(V,2)]
GL.VIEW(GL.numbering(.4)((V,[VV, EV])));


triangles = Lar.triangulate2d(V, EV)
GL.VIEW([ GL.GLGrid(V,triangles) ]);


# generate edges
ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
# remove duplicated edges from triangulation
nodupev = collect(Set(ev))
GL.VIEW(GL.numbering(.25)((V,[VV, nodupev, triangles])));
