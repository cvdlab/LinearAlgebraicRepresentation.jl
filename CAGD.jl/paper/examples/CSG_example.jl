# Rod Generation

npts = 16
V,_ = Lar.rod()([npts, 1])
EV = [
    [[2*i, (2*i+1)%2npts+1] for i = 1 : npts];          # horizontal upper edges
    [[2*i-1, (2*i)%2npts+1] for i = 1 : npts];          # horizontal lower edges
    [[2*i+1, 2*i+2] for i = 0 : npts-1];                # vertical edges
]
FE = [
    [i for i = 1 : npts],                                       # upper face
    [npts + i for i = 1 : npts],                                # lower face
    [[i, npts+i, 2*npts+i, 2*npts+i+1] for i = 1:npts-1]...,    # vertical faces
    [npts, 2*npts, 3*npts, 2*npts+1]
]
FV = Lar.cop2lar(Lar.lar2cop(FE) * Lar.lar2cop(EV))
CV = [collect(1 : size(V,2))]

# Three Rods Generation

cyl = Lar.Struct([ 
    Lar.r(0,0,0),
    Lar.s(1.0, 1.0, 1.5),
    Lar.t(0,0,-1.5),
    (V,CV,FV,EV)
])
cyl = Lar.struct2lar(cyl)
tris = Lar.Struct([
    cyl,
    Lar.Struct([Lar.r(pi/2,0,0), cyl ]) ,
    Lar.Struct([Lar.r(0,pi/2,0), cyl ])
])
V, CV, FV, EV = Lar.struct2lar(tris)
cyls = (V, CV, FV, EV)

# Cube Generation

V, (_, EV, FV, CV) = Lar.cuboid([1,1,1],true,[-1,-1,-1])
cube = (V, CV, FV, EV)

# Sphere Generation

V, FV = Lar.sphere(1.0)([15,25])
FV = sort(sort.(FV))
function getFaceEdges(fV::Array{Int,1})
    return [fV[1], fV[2]], [fV[1], fV[3]], [fV[2], fV[3]]
end
EVs = unique([(map(f -> getFaceEdges(f), FV)...)...])
CV = [collect(1:size(V, 2))]
sphere = (V, CV, FV, EVs)

# Object Creation

carry = Lar.Struct([
    cyls,
    Lar.s(1.3,1.3,1.3),
    cube,
    Lar.s(1.4,1.4,1.4),
    #Lar.s(1.5,1.5,1.5)
    sphere
])
V,CV,FV,EV = Lar.struct2lar(carry)

# Model Generation

model = CAGD.Model(V)
cFE = convert(Lar.ChainOp, Lar.coboundary_1(model.G, FV, EV))
CAGD.addModelCells!(model, 1, EV, signed = true)
CAGD.addModelCells!(model, 2, cFE)



A_model = CAGD.spatial_arrangement(model, outerCell = true)


V,CVs,FVs,EVs = Lar.pols2tria(A_model.G, A_model.T[1], A_model.T[2], A_model.T[3]);

if todisplay
    GL.VIEW(GL.GLExplode(V, FVs,  1.5,1.5,1.5,  99,1));
    GL.VIEW(GL.GLExplode(V, CVs[2:end],  5,5,5,  99,0.5));
end

##==============================================================================
##  Boolean Decomposition
##==============================================================================


boolean_matrix = Bool[
	1 0 0 0 0 0
	0 1 0 0 0 0
	0 1 0 0 0 0
	0 1 0 0 0 1
	0 0 0 0 0 1
	0 0 0 0 1 1

	0 1 0 0 1 1
	0 0 0 0 0 1
	0 1 0 0 0 1
	0 1 0 0 1 1
	0 0 0 1 1 1

	0 1 0 1 1 1
	0 1 0 1 1 1
	0 0 1 1 1 1
	0 0 1 1 1 1
	0 1 1 1 1 1

	0 0 1 0 1 1
	0 1 1 0 1 1
	0 1 1 0 1 1
	0 0 1 1 1 1
	0 0 0 1 1 1

	0 1 0 1 1 1
	0 1 0 1 1 1
	0 0 1 1 1 1
	0 0 1 0 1 1
	0 1 1 0 1 1

	0 1 1 0 1 1
	0 0 1 0 0 0
	0 0 1 0 0 0
	0 0 1 0 0 1
	0 0 0 0 0 1

	0 0 1 0 0 1
	0 0 0 0 0 1
	0 0 0 1 0 0
	0 0 0 1 0 0
	0 0 0 1 0 1

	0 0 0 0 0 1
	0 0 0 1 0 1
	0 0 0 0 0 1
	0 0 0 0 1 0
]

bm = sparse(boolean_matrix);
external = bm[:,1];
rod1     = bm[:,2];
rod2     = bm[:,3];
rod3     = bm[:,4];
cube     = bm[:,5];
sphere   = bm[:,6];

shells = Array{SparseVector{Bool,Int64}, 1}();

push!(shells, .|(sphere, cube, rod1, rod2, rod3));
push!(shells, .|(rod1, rod2));
push!(shells, rod1);

push!(shells, .&(sphere, cube));
push!(shells, .&(sphere, .!(cube, rod1, rod2, rod3)));
push!(shells, .&(sphere, cube, .!rod1, .!rod2, .!rod3));

push!(shells, .|(shells[6], shells[3]));
push!(shells, .|(rod1, rod2, rod3));
push!(shells, .|(sphere, cube, rod1, rod2, rod3));

if todisplay  for i = 1 : 9
    GL.VIEW(GL.GLExplode(V,CVs[shells[i].nzind], 1, 1, 1, 99, 0.4));
    GL.VIEW(GL.GLExplode(V,CVs[shells[i].nzind], 4, 4, 4, 99, 0.4));
end  end