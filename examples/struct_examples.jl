using LARLIB
using LARVIEW

square=([[0; 0] [0; 1] [1; 0] [1; 1]], [[1, 2, 3, 4]], 
[[1,2], [1,3], [2,4], [3,4]])
V,FV,EV = square
model = V,([[1],[2],[3],[4]],EV,FV)
table=LARLIB.apply(LARLIB.t(-0.5,-0.5), square)
chair=LARLIB.Struct([LARLIB.t(0.75,0),LARLIB.s(0.35,0.35),table])
structo=LARLIB.Struct([LARLIB.t(2,1),table,repeat([LARLIB.r(pi/2),chair],outer=4)...])
structo1=LARLIB.Struct(repeat([structo,LARLIB.t(0,2.5)],outer=10));
structo2=LARLIB.Struct(repeat([structo1,LARLIB.t(3,0)],outer=10));
scene=LARLIB.evalStruct(structo2);
LARVIEW.view(scene)
W,FW,EW = LARLIB.struct2lar(structo2);
LARVIEW.view(LARVIEW.lar2hpc(W,EW))
assembly = LARLIB.Struct([LARLIB.sphere()(), LARLIB.t(3,0,-1), LARLIB.cylinder()()])
hpc = LARVIEW.mkpol(LARLIB.struct2lar(assembly)...)
LARVIEW.view(hpc)
LARVIEW.view(LARLIB.struct2lar(assembly))

