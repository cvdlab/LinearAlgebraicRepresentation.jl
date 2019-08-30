using LinearAlgebraicRepresentation
L = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL
# using Revise

# square=([[0; 0] [0; 1] [1; 0] [1; 1]], [[1, 2, 3, 4]],
# [[1,2], [1,3], [2,4], [3,4]])
# V,FV,EV = square
# model = V,([[1],[2],[3],[4]],EV,FV)
# table=L.apply(L.t(-0.5,-0.5), square)
# chair=L.Struct([L.t(0.75,0),L.s(0.35,0.35),table])
# structo=L.Struct([L.t(2,1),table,repeat([L.r(pi/2),chair],outer=4)...])
# structo1=L.Struct(repeat([structo,L.t(0,2.5)],outer=10));
# structo2=L.Struct(repeat([structo1,L.t(3,0)],outer=10));
#
# V,FV,EV = L.struct2lar(structo2);
# GL.VIEW([ GL.GLHulls(V,FV, GL.COLORS[1]) ]);
# #GL.VIEW([ GL.GLHulls(V,EV, GL.COLORS[1]) ]);  # BUG with colors?
#
# assembly = (L.struct2lar ∘ L.Struct)([L.sphere()(), L.t(3,0,-1), L.cylinder()()])
# GL.VIEW([ GL.GLHulls(L.cylinder()()...,GL.Point4d(1,1,1,0.2)) ]);
# GL.VIEW([ GL.GLHulls(assembly...,GL.Point4d(1,1,1,0.2)) ]);

cube = L.apply( L.t(-.5,-.5,0), L.cuboid([1,1,1]))
tableTop = L.Struct([ L.t(0,0,.85), L.s(1,1,.05), cube ])
tableLeg = L.Struct([ L.t(-.475,-.475,0), L.s(.1,.1,.89), cube ])
tablelegs = L.Struct( repeat([ tableLeg, L.r(0,0,pi/2) ],outer=4) )
table = L.Struct([ tableTop, tablelegs ])
table = L.struct2lar(table)
GL.VIEW([ GL.GLPol(table..., GL.Point4d(0,1,0,0.1)) ]);

cylndr = L.rod(.06, .5, 2*pi)([8,1])
chairTop = L.Struct([ L.t(0,0,0.5), L.s(0.5,0.5,0.04), cube ])
chairLeg = L.Struct([ L.t(-.22,-.22,0), L.s(.5,.5,1), L.r(0,0,pi/8), cylndr ])
chairlegs = L.Struct( repeat([ chairLeg, L.r(0,0,pi/2) ],outer=4) );
chair = L.Struct([ chairTop, chairlegs ]);
chair = L.struct2lar(chair)
GL.VIEW([ GL.GLPol(chair..., GL.Point4d(1,1,1,0.2)) ]);

theChair = L.Struct([ L.t(-.8,0,0), chair ])
fourChairs = L.Struct( repeat([L.r(0,0,pi/2), theChair],outer=4) );
fourSit = L.Struct([fourChairs,table]);
fourSit = L.struct2lar(fourSit)
GL.VIEW([ GL.GLPol(fourSit..., GL.Point4d(1,1,1,0.2)) ]);

singleRow=L.Struct(repeat([fourSit,L.t(0,2.5,0)],outer=10));
singleRow = L.struct2lar(singleRow)
GL.VIEW([ GL.GLPol(singleRow..., GL.Point4d(1,1,1,0.2)) ]);

refectory=L.Struct(repeat([singleRow,L.t(3,0,0)],outer=10));
refectory = L.struct2lar(refectory)
GL.VIEW([ GL.GLPol(refectory..., GL.Point4d(1,1,1,0.2)) ]);
