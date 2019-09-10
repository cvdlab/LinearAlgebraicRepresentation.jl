using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL
using Base.union

# 3D Boolean example generation
#-------------------------------------------------------------------------------
n,m,p = 1,1,1
V,(VV,EV,FV,CV) = Lar.cuboidGrid([n,m,p],true)
cube = V,FV,EV

V,FV = Lar.sphere()()
EV = Lar.simplexFacets(FV)
sphere = V,FV,EV

# three cubes in "assembly"
assembly = Lar.Struct([ cube,
    Lar.t(.3,.4,.25), Lar.r(pi/5,0,0), Lar.r(0,0,pi/12), cube,
    Lar.t(-.2,.4,-.2), Lar.r(0,pi/5,0), Lar.r(0,pi/12,0), cube ])

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
AxorB = AorB .& (.! AandB)
AorBorC = A .| B .| C
AorBorC = .|(A,B,C)
AandBandC = A .& B .& C
AandBandC = .&(A,B,C)
AminusBminusC = .&(A, .!(B .| C)) # A - B - C

union = .|(A,B,C)
union = Matrix(copCF)' * Int.(union) # coord vector of Faces
intersection = Matrix(copCF)' * Int.(AandBandC) # coord vector of Faces
difference = Matrix(copCF)' * Int.(AminusBminusC) # coord vector of Faces

V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF) # whole assembly
Fs = union
V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF, Fs) # part of assembly





# EV = Lar.cop2lar(copEV)
# EVor = [ev for (k,ev) in enumerate(EV) if abs(union[k])==1 ]
# EVand = [ev for (k,ev) in enumerate(EV) if abs(intersection[k])==1 ]
# EVxor = [ev for (k,ev) in enumerate(EV) if abs(xor[k])==1 ]




GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
GL.VIEW(GL.GLExplode(V,EVs,1.,1.,1.,99,1));
meshes = GL.GLExplode(V,CVs[2:end],1.5,1.5,1.5,99,1);
GL.VIEW( push!( meshes, GL.GLFrame) );

# GL.VIEW(GL.GLExplode(V,[EVor],1.,1.,1.,99,1));


GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],0.5), GL.GLFrame ]);
