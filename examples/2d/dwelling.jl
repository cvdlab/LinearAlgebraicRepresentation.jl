using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

using Revise, OhMyREPL

wall = Lar.svg2lar("test/svg/assembly/wall.svg", flag=false)
openings = Lar.svg2lar("test/svg/assembly/openings.svg", flag=false)
rectangle = Lar.svg2lar("test/svg/assembly/rectangle.svg", flag=false)
box = Lar.svg2lar("test/svg/assembly/box.svg", flag=false)

assembly = Lar.Struct([wall , openings, rectangle, box])
V,EV = Lar.struct2lar(assembly)
GL.VIEW([ GL.GLFrame2, GL.GLGrid(V,EV, GL.COLORS[1],1) ]);
W, copEV, copFE, boolmatrix = Lar.bool2d(assembly)
@show Matrix(boolmatrix)

innerpoints = Lar.internalpoints2d(W,copEV,copFE[1:end,:])
points = convert(Array{Float64,2}, hcat(innerpoints...)')
GL.VIEW([ GL.GLFrame2, GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLPoints(points) ]);


A = boolmatrix[:,1]
B = boolmatrix[:,2]
C = boolmatrix[:,3]
D = boolmatrix[:,4]
expr = ( D .& .!(.&(A, .!B, .!C)) )# A - B - C
formula = Matrix(copFE)' * Int.(expr) # coord vector of Faces
V = convert(Lar.Points,W')
EV = Lar.cop2lar(copEV)
EVmin = [ev for (k,ev) in enumerate(EV) if abs(formula[k])==1 ]
V = Lar.normalize(V)
V = Lar.normalize(V)
GL.VIEW([ GL.GLFrame2, GL.GLGrid(V,EVmin, GL.COLORS[1],1) ]);

VV = [[k] for k in 1:size(V,2)];
model = (V,[VV, EV])
GL.VIEW( GL.numbering(.1)( model,GL.COLORS[1],0.1 ) );
