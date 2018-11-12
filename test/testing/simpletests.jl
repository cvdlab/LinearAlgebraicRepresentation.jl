using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

#translate = Lar.t(0,0,0) #OK
#translate = Lar.t(0.25,0,0) #OK
#translate = Lar.t(0.25,0,0.25) #OK
#translate = Lar.t(0.25,0,-0.25) #OK
#translate = Lar.t(.25,.25,1.0) #OK
#translate = Lar.t(.25,1,-.25) #OK
#translate = Lar.t(.25,0.9,-.25) #OK
translate = Lar.t(.25,.25,-.25) #K0

V1,(_,EV1,FV1,CV1) = Lar.cuboid([1,1,1],true)
V2,(_,EV2,FV2,CV2) = Lar.cuboid([.25,.5,1.5],true)
assembly = Lar.Struct([(V1,FV1,EV1),translate,(V2,FV2,EV2)])
W,FW,EW = Lar.struct2lar(assembly)
Plasm.view(W,EW)
V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW )
Plasm.view(V,EV)
