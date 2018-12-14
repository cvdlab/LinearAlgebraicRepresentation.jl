using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using Plasm

transf = Lar.t(0,0,0) #OK
transf = Lar.t(-0.125,-0.125,-0.125) #0K
transf = Lar.t(0.25,0,0) #OK
transf = Lar.t(0.25,0,0.25) #OK
transf = Lar.t(0.25,0,-0.25) #OK
transf = Lar.t(.25,.25,1.0) #OK
transf = Lar.t(.25,1.,-.25) #OK
transf = Lar.t(.25,0.9,-.25) #	KO !!!!
transf = Lar.t(.25,.25,-.25) #0K
transf = Lar.t(.25,.25,0) #0K
transf = Lar.t(0.75,0.75,-.25) #0K
transf = Lar.t(0.875,0.875,-.25) #0K
transf = Lar.t(0.5,0.5,0)*Lar.r(0,0,pi/3)*Lar.t(-0.125,-0.25,0) #OK

V1,(_,EV1,FV1,CV1) = Lar.cuboid([1,1,1],true)
V2,(_,EV2,FV2,CV2) = Lar.cuboid([.25,.5,1.5],true)
assembly = Lar.Struct([(V1,FV1,EV1),transf,(V2,FV2,EV2)])
W,FW,EW = Lar.struct2lar(assembly)
Plasm.view(W,EW)
V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW )
Plasm.view(Plasm.numbering(.25)((V,[[[k] for k=1:size(V,2)],EV,FV]))) 






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
Lar.t(.25,.25,-.25),
Lar.t(1.125,1.125,-.25),
Lar.t(0.75,0.75,-.25),
Lar.t(0.875,0.875,-.25),
Lar.t(0.5,0.5,0)*Lar.r(0,0,pi/3)*Lar.t(-0.125,-0.25,0) #KO
]
 
for transf in maps
	V1,(_,EV1,FV1,CV1) = Lar.cuboid([1,1,1],true)
	V2,(_,EV2,FV2,CV2) = Lar.cuboid([.25,.5,1.5],true)
	assembly = Lar.Struct([(V1,FV1,EV1),transf,(V2,FV2,EV2)]);
	W,FW,EW = Lar.struct2lar(assembly);
	#Plasm.view(W,EW)
	V,(EV,FV,CV),(copEV,copFE,copCF) = Lar.chaincomplex( W,FW,EW );
	@show CV
	Plasm.view(Plasm.numbering(.25)((V,[[[k] for k=1:size(V,2)],EV,FV]))) 
end

