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
Plasm.view(Plasm.numbering(0.35)((V,[VV, EV])))

triangles = Lar.triangulate2d(V, EV)
Plasm.view(V,triangles)

ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
# remove duplicated edges from triangulation
ev_nodups = collect(Set(ev))
Plasm.view(Plasm.numbering(0.25)((V,[VV, ev_nodups, triangles])))
