using LinearAlgebraicRepresentation, ViewerGL, SparseArrays
Lar = LinearAlgebraicRepresentation; GL = ViewerGL

wall = Lar.svg2lar("test/svg/assembly/wall.svg", flag=false)
openings = Lar.svg2lar("test/svg/assembly/openings.svg", flag=false)
rectangle = Lar.svg2lar("test/svg/assembly/rectangle.svg", flag=false)
box = Lar.svg2lar("test/svg/assembly/box.svg", flag=false)

assembly = Lar.Struct([wall, openings, rectangle])
V,EV = Lar.boolops(assembly, :|)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = Lar.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

assembly = Lar.Struct([wall, openings, rectangle])
V,EV = Lar.boolops(assembly, :-)
diff = Lar.Struct([ (V,EV) ])
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = Lar.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

assembly = Lar.Struct([wall, openings, box])
V,EV = Lar.boolops(assembly, :&)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = Lar.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

inner = Lar.Struct([ box,diff ])
V,EV = Lar.boolops( inner, :-)
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);
V,FVs,EVs = Lar.arrange2D(V,EV)
GL.VIEW(GL.GLExplode(V,FVs,1.,1.,1.,99,1));

W, copEV, copFE, boolmatrix = Lar.bool2d(inner)
FVs = Lar.triangulate2D(W, [copEV, copFE])
V = convert(Lar.Points, W')
GL.VIEW(GL.GLExplode(V,FVs,1.5,1.5,1.5,99,1));
