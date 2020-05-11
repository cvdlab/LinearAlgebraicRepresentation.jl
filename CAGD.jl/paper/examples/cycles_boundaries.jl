npts = 64;

## METHOD 1

function cylinder(;r=1, h=1, n=8, k=1)
	V,FV = Lar.ring(0,r)([n,r])
	W,CW = Lar.extrudeSimplicial((V,FV), [h for i=1:k])
	V = W .- [0,0,h*k/2]
	FW = Lar.simplexFacets(CW)
	FV = Lar.simplexBoundary_3(CW,FW)
	EV = [[u,v] for (u,v) in sort(collect(Set(map(sort,cat(
					[[[u,v],[v,w],[w,u]] for (u,v,w) in FV]))))) if u!=v]
	return V,FV,EV
end

tube = cylinder(n=npts,h=2,k=2);

triple = Lar.Struct([ tube, Lar.r(pi/2,0,0), tube, Lar.r(0,pi/2,0), tube ]);
model = Lar.struct2lar(triple)
V,FV,EV = model
GL.VIEW([
	GL.GLFrame,
	GL.GLGrid( model[1:2]..., GL.COLORS[1], 0.5)
]);


model = CAGD.Model(V);
cFE = convert(Lar.ChainOp, Lar.coboundary_1(model.G, FV, EV));
CAGD.addModelCells!(model, 1, EV, signed = true);
CAGD.addModelCells!(model, 2, cFE);

A_model = CAGD.spatial_arrangement(model);

if todisplay
    V, CVs, FVs, EVs = Lar.pols2tria(A_model.G, A_model.T[1], A_model.T[2], A_model.T[3]);

    GL.VIEW(GL.GLExplode(V, FVs,  1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, EVs,  1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, CVs[1:end],  5,5,5,  99,0.5));
end;




## METHOD 2

# Build a cylinder

V,_ = Lar.rod()([npts, 1]);
EV = [
    [[2*i, (2*i+1)%2npts+1] for i = 1 : npts];          # horizontal upper edges
    [[2*i-1, (2*i)%2npts+1] for i = 1 : npts];          # horizontal lower edges
    [[2*i+1, 2*i+2] for i = 0 : npts-1];                # vertical edges
];
FE = [
    [i for i = 1 : npts],                                       # upper face
    [npts + i for i = 1 : npts],                                # lower face
    [[i, npts+i, 2*npts+i, 2*npts+i+1] for i = 1:npts-1]...,    # vertical faces
    [npts, 2*npts, 3*npts, 2*npts+1]
];

FV = Lar.cop2lar(Lar.lar2cop(FE) * Lar.lar2cop(EV));


# Build three cylinders

cyl = Lar.Struct([
    Lar.r(0,0,0),
    Lar.s(1.0, 1.0, 1.5),
    Lar.t(0,0,-1.5),
    (V,FV,EV)
]);
cyl = Lar.struct2lar(cyl);
tris = Lar.Struct([
    cyl,
    Lar.Struct([Lar.r(pi/2,0,0), cyl ]) ,
    Lar.Struct([Lar.r(0,pi/2,0), cyl ])
]);
V, FV, EV = Lar.struct2lar(tris);

if todisplay
    GL.VIEW([
        GL.GLFrame
        GL.GLGrid(V,EV,GL.COLORS[1],1)
    ]);
end;

# Model Generation

model = CAGD.Model(V);
cFE = convert(Lar.ChainOp, Lar.coboundary_1(model.G, FV, EV));
CAGD.addModelCells!(model, 1, EV, signed = true);
CAGD.addModelCells!(model, 2, cFE);

A_model = CAGD.spatial_arrangement(model);

if todisplay
    V, CVs, FVs, EVs = Lar.pols2tria(A_model.G, A_model.T[1], A_model.T[2], A_model.T[3])

    GL.VIEW(GL.GLExplode(V, FVs,  1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, EVs,  1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, CVs[1:end],  5,5,5,  99,0.5));
end;
