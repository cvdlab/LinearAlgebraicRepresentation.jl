# Squares building

W,(WW,EW,FW) = Lar.cuboid([1,1], true)
square = (W,EW)

squares = Lar.Struct([  square,  Lar.t(-1, -1),  Lar.s(3, 3),  square  ])

V,EV = Lar.struct2lar(squares)

model = CAGD.Model(V);
CAGD.addModelCells!(model, 1, EV);

# Lines adding
lines = CAGD.Model([
    -2.0  1.0  0.0  3.0
     0.0  3.0  3.0  0.0
]);
CAGD.addModelCells!(lines, 1, [[1, 2], [3, 4]]);


# Merge the models
CAGD.uniteModels!(model, lines)


# Planar Arrangment

A_model = CAGD.planar_arrangement(model, sparse(Int8[0,0,0,0,1,1,1,1,0,0]));

A_model.G
Matrix(A_model.T[1])
Matrix(A_model.T[2])

# Display Models

if todisplay
    VV = [[k] for k=1:size(A_model.G ,2)];
    triangulated_faces = CAGD.triangulate2D(A_model);
    FVs = convert(Array{Lar.Cells}, triangulated_faces);

    GL.VIEW(GL.GLExplode(A_model.G, FVs, 1.0, 1.0, 1.0, 99, 1));
    GL.VIEW(GL.GLExplode(A_model.G, FVs, 1.2, 1.2, 1.2, 99, 1));
end
