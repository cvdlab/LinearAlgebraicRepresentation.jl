# Generate a 100 grid of squares

V, EV = Lar.randomcuboids(100, .4);
V = GL.normalize2(V, flag = true);

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
end;

# Generate 50 S1

store = []
mycircle(r,n) = Lar.circle(r)([n])
for k=1:50
    r = rand()
    n = 8
    while true
        n = abs(rand(Int8)+1)
        if n>2 break end
    end
    scale = Lar.s(0.25,0.25)
    transl = Lar.t(rand(2)...)
    s = Lar.Struct([ transl, scale, mycircle(r,n) ])
    push!(store, Lar.struct2lar(s))
end

V,EV = Lar.struct2lar(Lar.Struct(store))
V = GL.normalize2(V)

if todisplay  GL.VIEW([ GL.GLLines(V,EV) ])  end

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
end;
