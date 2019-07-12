using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL
using Triangle, LinearAlgebra

# input of primitive shapes
V,(VV,EV,FV) = Lar.simplex(2, true)
W,(WW,EW,FW) = Lar.cuboid([1,1], true)
triangle = (V,EV)
square = (W,EW)
# hierarchical assembly
model = Lar.Struct([
            Lar.Struct([
                square,
                Lar.t(.175,.175), Lar.s(.5,.5),
                triangle
            ]),
            Lar.t(.825,.825), Lar.s(-.5,-.5),
            triangle
        ])
# visualization of generated wire-frame model
V,EV = Lar.struct2lar(model)
GL.VIEW(GL.numbering(.4)((V,[VV, EV])));

triangles = Lar.triangulate2d(V, EV)
GL.VIEW([ GL.GLGrid(V,triangles) ]);

# generate edges
ev = map(sort,cat([[[u,v], [v,w], [w,u]] for (u,v,w) in triangles]))
# remove duplicated edges from triangulation
nodupev = collect(Set(ev))

GL.VIEW(GL.numbering(.3)((V,[VV, nodupev])));
GL.VIEW([ GL.GLGrid(V,triangles) ]);
GL.VIEW(GL.numbering(.25)((V,[VV, nodupev, triangles]),GL.COLORS[1],0.2));
