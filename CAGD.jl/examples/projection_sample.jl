using CAGD;

model = CAGD.Model([
    2.0 1.0 3.0 2.0 1.0 -3.
    2.0 4.0 5.0 1.0 9.0 1.0
    1.0 4.0 1.0 4.0 3.0 4.0
]);

CAGD.addModelCells!(model, 1, [[1,2],[2,3],[1,3],[1,4],[2,4],[1,5],[2,5],[1,6],[2,6]]);
CAGD.addModelCells!(model, 2, [[1,2,3],[1,4,5],[1,6,7],[1,8,9]]);

#test = eval_ord_faces(model, 1)
τ = 1;
faces = model.T[2][:, τ].nzind
edge = model.T[1][τ, :].nzind
FVs = [setdiff((model.T[1]'*model.T[2][f, :]).nzind, edge) for f in faces]
Gface, Gfaceidx = CAGD.getModelCellGeometry(model, 2, faces[1], true)
M = CAGD.build_projection_matrix(Gface, y0 = true)
PG = (M * [model.G; ones(1, size(model, 0, 2))])[2:3, :]
angles = [[mod(atan(PG[2, v], PG[1, v]), 2π) for v in fv] for fv in FVs]
angles = [[α for α in a if !(α ≈ 0)] for a in angles]
angles = [[faces[i], angles[i][1]] for i = 1 : length(faces)]
sort!(angles, lt = (a, b) -> a[2] < b[2])
map(a->Int(a[1]), angles)