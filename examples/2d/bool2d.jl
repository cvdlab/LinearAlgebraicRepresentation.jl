using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

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



n,m = 1,1
V,(VV,EV,FV) = Lar.cuboidGrid([n,m],true)
square = V,EV
q = sqrt(2)
V, (VV,EV,FV) = Lar.cuboid([-q,-q], true, [q,q])
assembly = Lar.Struct([ #(V,FV,EV),
    Lar.Struct([ Lar.t(-sqrt(2)/4, -sqrt(2)/2 ), Lar.r(pi/4), square ]),
    Lar.Struct([ Lar.t( sqrt(2)/4, -sqrt(2)/2 ), Lar.r(pi/4), square ])
])
V,EV = Lar.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);

# VV = [[k] for k in 1:size(V,2)];
# model = (V,[VV, EV])
# GL.VIEW( GL.numbering(.75)( model,GL.COLORS[1],0.1 ) );

# booleanmatrix = Lar.booleanops2d((V,EV))
booleanmatrix = Lar.booleanops2d(assembly)
@show booleanmatrix

# using BitBasis
#
# julia> a = BitArray([0,1,1,1,0,1])
# 6-element BitArray{1}:
# false
#  true
#  true
#  true
# false
#  true
#
# julia> b = map(! , BitArray([0,1,1,1,0,1]))
# 6-element BitArray{1}:
#  true
# false
# false
# false
#  true
# false
#
# julia> a .& b
# 6-element BitArray{1}:
# false
# false
# false
# false
# false
# false
#
# julia> a .| b
# 6-element BitArray{1}:
# true
# true
# true
# true
# true
# true
