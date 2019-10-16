using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
import Base.union

# 3D Boolean example generation
#-------------------------------------------------------------------------------
n,m,p = 1,1,1
V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
cube = V,FV,EV

# three cubes in "assembly"
assembly = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube
])

V,FV,EV = Lar.struct2lar(assembly)
GL.VIEW([ GL.GLGrid(V,FV), GL.GLFrame ]);

W, (copEV, copFE, copCF), boolmatrix = Lar.bool3d(assembly)
Matrix(boolmatrix)
#three-chains = [ for k = 1:3]

A = boolmatrix[:,2]
B = boolmatrix[:,3]
C = boolmatrix[:,4]
AorB = A .| B
AandB = A .& B
AxorB = AorB .& (.!AandB) # = A .⊻ B
AorBorC = A .| B .| C
AorBorC = .|(A, B, C)
AandBandC = A .& B .& C
AandBandC = .&(A, B, C)
AminBminC = .&(A, .!B, .!C) # A - B - C

unione = Matrix(copCF)' * Int.(AorBorC) # coord vector of Faces
intersection = Matrix(copCF)' * Int.(AandBandC) # coord vector of Faces
difference = Matrix(copCF)' * Int.(AminBminC) # coord vector of Faces

V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF) # whole assembly
Fs = difference
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF, Fs) # part of assembly


# EV = Lar.cop2lar(copEV)
# EVor = [ev for (k,ev) in enumerate(EV) if abs(unione[k])==1 ]
# EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
# EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]


GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,1,1));
# meshes = GL.GLExplode(V,CVs,5,5,5,99,1);
# GL.VIEW( push!( meshes, GL.GLFrame) );
