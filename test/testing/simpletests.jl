using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

#translate = Lar.t(0,0,0) #OK
#translate = Lar.t(-0.125,-0.125,-0.125) #K0
#translate = Lar.t(0.25,0,0) #OK
#translate = Lar.t(0.25,0,0.25) #OK
#translate = Lar.t(0.25,0,-0.25) #OK
#translate = Lar.t(.25,.25,1.0) #OK
#translate = Lar.t(.25,1.,-.25) #OK
translate = Lar.t(.25,0.9,-.25) #OK
#translate = Lar.t(.25,.25,-.25) #0K

V1,(_,EV1,FV1,CV1) = Lar.cuboid([1,1,1],true)
V2,(_,EV2,FV2,CV2) = Lar.cuboid([.25,.5,1.5],true)
assembly = Lar.Struct([(V1,FV1,EV1),translate,(V2,FV2,EV2)])
W,FW,EW = Lar.struct2lar(assembly)
Plasm.view(W,EW)
V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW )
Plasm.view(V,EV)







using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

maps = [
Lar.r(0,0,pi/4), 
Lar.t(-0.125,-0.125,-0.125), 
Lar.t(0.25,0,0), 
Lar.t(0.25,0,0.25), 
Lar.t(0.25,0,-0.25), 
Lar.t(.25,.25,1.0), 
Lar.t(.25,1.,-.25), 
Lar.t(.25,0.9,-.25), 
Lar.t(.25,.25,-.25)
]
 
for translate in maps
	V1,(_,EV1,FV1,CV1) = Lar.cuboid([1,1,1],true)
	V2,(_,EV2,FV2,CV2) = Lar.cuboid([.25,.5,1.5],true)
	assembly = Lar.Struct([(V1,FV1,EV1),translate,(V2,FV2,EV2)]);#,Lar.r(0,0,-pi/6),(V2,FV2,EV2)]);
	W,FW,EW = Lar.struct2lar(assembly);
	Plasm.view(Plasm.numbering(.25)((W,[[[k] for k=1:size(W,2)],EW,FW]))) 
	#Plasm.view(W,EW)
	V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW );
	Plasm.view(V,EV)
end


Plasm.view(Plasm.numbering(.25)((W,[[[k] for k=1:size(W,2)],EW,FW]))) 