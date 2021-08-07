wall = Lar.svg2lar("./test/svg/assembly/wall.svg", flag=false);
openings = Lar.svg2lar("./test/svg/assembly/svg/openings.svg", flag=false);
rectangle = Lar.svg2lar("./test/svg/assembly/rectangle.svg", flag=false);
box = Lar.svg2lar("./test/svg/assembly/box.svg", flag=false);

assembly = Lar.Struct([wall, openings, rectangle, box]);
V,EV = Lar.struct2lar(assembly);
if todisplay  GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],1), GL.GLFrame2 ]);  end

model = CAGD.Model(V);
CAGD.addModelCells!(model, 1, EV);

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
    GL.VIEW(GL.GLExplode(A_model.G, FVs, 4.0, 4.0, 4.0, 99, 1));
end

innerpoints = Lar.internalpoints2d(A_model.G', A_model.T[1], A_model.T[2][1:end, :])
# associate internal points to 2-cells
listOfModels = Lar.evalStruct(assembly)
inputfacenumbers = [length(listOfModels[k][2]) for k=1:length(listOfModels)];
# test input data for containment of reference points
boolmatrix = sparse(BitArray(undef, length(innerpoints), length(listOfModels)));
containmenttest = Lar.testinternalpoint2d(listOfModels);
for (k,point) in enumerate(innerpoints) # k runs on columns
    cells = containmenttest(point); # contents of columns
    #println(k," ",faces)
    for l in cells
        boolmatrix[k,l] = 1;
    end
end

boolmatrix = sparse(Bool[
    true  false false true
    true  false true  true
    false false false true
    true  true  false true
    false true  false true
    false false false true
    true  true  false true
    false false false true
    false true  false false
    false false true  true
    true  false false true
    false true  false true
    true  false false true
    true  false false true
    false true  false true
    false false  true  true
    false false true  false
]);

Matrix(boolmatrix)

wall = boolmatrix[:, 1]
Lshape = boolmatrix[:, 2]
rect = boolmatrix[:, 3]
box = boolmatrix[:, 4]

GL.VIEW(GL.GLExplode(A_model.G, FVs[wall.nzind],   1.0, 1.0, 1.0, 99, 1));
GL.VIEW(GL.GLExplode(A_model.G, FVs[Lshape.nzind], 1.0, 1.0, 1.0, 99, 1));
GL.VIEW(GL.GLExplode(A_model.G, FVs[rect.nzind],   1.0, 1.0, 1.0, 99, 1));
GL.VIEW(GL.GLExplode(A_model.G, FVs[box.nzind],    1.0, 1.0, 1.0, 99, 1));

shells = Array{SparseVector{Bool,Int64},1}()
push!(shells, .|(wall, Lshape, rect));
push!(shells, .&(wall, .!Lshape, .!rect));
push!(shells, .&(box, .!shells[2]))

if todisplay
    for i = 1 : length(shells)
    GL.VIEW(GL.GLExplode(A_model.G, FVs[shells[i].nzind], 1.0, 1.0, 1.0, 1, 1));
    GL.VIEW(GL.GLExplode(A_model.G, FVs[shells[i].nzind], 1.0, 1.0, 1.0, 99, 1));
    end
end
