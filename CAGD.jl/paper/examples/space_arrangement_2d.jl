W,(WW,EW,FW) = Lar.cuboid([1,1], true)
square = (W,EW)

square = Lar.Struct([  Lar.r(pi/4),  square  ])
squares = Lar.Struct([  square,  Lar.t(√2/2, 0),  square  ])

V,EV = Lar.struct2lar(squares)
VV = [[k] for k=1:size(V,2)]
GL.VIEW(GL.numbering(.4)((V,[VV, EV])));

model = CAGD.Model(V);
CAGD.addModelCells!(model, 1, EV);

# Planar Arrangment

A_model = CAGD.planar_arrangement(model);

A_model.G
Matrix(A_model.T[1])
Matrix(A_model.T[2])

# Display Models

if todisplay
    VV = [[k] for k=1:size(A_model.G, 2)];
    triangulated_faces = CAGD.triangulate2D(A_model);
    FVs = convert(Array{Lar.Cells}, triangulated_faces);

    GL.VIEW(GL.GLExplode(A_model.G, FVs, 1.0, 1.0, 1.0, 99, 1));
    GL.VIEW(GL.GLExplode(A_model.G, FVs, 1.2, 1.2, 1.2, 99, 1));
end

innerpoints = Lar.internalpoints2d(A_model.G', A_model.T[1], A_model.T[2][1:end, :])
# associate internal points to 2-cells
listOfModels = Lar.evalStruct(squares)
inputfacenumbers = [length(listOfModels[k][2]) for k=1:length(listOfModels)]
# test input data for containment of reference points
boolmatrix = BitArray(undef, length(innerpoints), length(listOfModels))
containmenttest = Lar.testinternalpoint2d(listOfModels)
for (k,point) in enumerate(innerpoints) # k runs on columns
    cells = containmenttest(point) # contents of columns
    #println(k," ",faces)
    for l in cells
        boolmatrix[k,l] = 1
    end
end

sq1 = boolmatrix[:, 1]
sq2 = boolmatrix[:, 2]

sqXor = sq1 .⊻ sq2

if todisplay
    GL.VIEW(GL.GLExplode(A_model.G, FVs[sparse(sqXor).nzind], 1.0, 1.0, 1.0, 99, 1));
end

