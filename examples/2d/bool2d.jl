using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

using Revise, OhMyREPL


# 3D Boolean example generation
#-------------------------------------------------------------------------------
n,m = 1,1
V,(VV,EV,FV) = Lar.cuboidGrid([n,m],true)
square = V,FV,EV
q = sqrt(2)
V, (VV,EV,FV) = Lar.cuboid([-q,-q], true, [q,q])

twosquares = Lar.Struct([ #(V,FV,EV),
    Lar.Struct([ Lar.t(-sqrt(2)/4, -sqrt(2)/2 ), Lar.r(pi/4), square ]),
    Lar.Struct([ Lar.t( sqrt(2)/4, -sqrt(2)/2 ), Lar.r(pi/4), square ])
    ])

V,FV,EV = Lar.struct2lar(twosquares)
GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame2 ]);
W,EW = Lar.fragmentlines((V,EV))
V,FVs,EVs = Lar.arrange2D(W,EW)
GL.VIEW(GL.GLExplode(V,FVs[1:3],1.,1.,1.,99,1));

circle = Lar.circle(1.06066/2)()

n,m = 1,1
V,(VV,EV,FV) = Lar.cuboidGrid([n,m],true)
square = V,EV
q = sqrt(2)
assembly = Lar.Struct([ #(V,FV,EV),
    Lar.Struct([ Lar.t(-q/4, 0 ), Lar.r(pi/4), circle ]),
    Lar.Struct([ Lar.t( q/4, -q/2 ), Lar.r(pi/4), square ])
])
V,EV = Lar.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

# boolmatrix = Lar.booleanops2d((V,EV))
W, copEV, copFE, boolmatrix = Lar.bool2d(assembly)
@show Matrix(boolmatrix)
@show Matrix(copFE)
A = boolmatrix[:,1]
B = boolmatrix[:,2]
AorB = A .| B
AandB = A .& B
AxorB = AorB .& (.! AandB)

union = Matrix(copFE)' * Int.(AorB)
intersection = Matrix(copFE)' * Int.(AandB)
xor = Matrix(copFE)' * AxorB

V = convert(Lar.Points,W')
EV = Lar.cop2lar(copEV)
EVor = [ev for (k,ev) in enumerate(EV) if abs(union[k])==1 ]
EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]
GL.VIEW([ GL.GLGrid(V,EVor, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVand, GL.COLORS[1],1), GL.GLFrame2 ]);
GL.VIEW([ GL.GLGrid(V,EVxor, GL.COLORS[1],1), GL.GLFrame2 ]);

VV = [[k] for k in 1:size(V,2)];
model = (V,[VV, EV])
GL.VIEW( GL.numbering(.5)( model,GL.COLORS[1],0.1 ) );
