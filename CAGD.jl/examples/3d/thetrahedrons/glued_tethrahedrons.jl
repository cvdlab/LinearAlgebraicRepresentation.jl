V = [
    0.0 +1  0  0  0
    0.0  0  0 +1 -1
    0.0  0 +1  0  0
];

W = [
     0.0  0.0  0.0
    +1.0  0.0  0.0
     0.0  0.0 +1.0
     0.0 +1.0  0.0
     0.0 -1.0  0.0
];

EV = SparseArrays.sparse(Int8.([
    -1 +1  0  0  0
    -1  0 +1  0  0
    -1  0  0 +1  0
    -1  0  0  0 +1
     0 -1 +1  0  0
     0 -1  0 +1  0
     0 -1  0  0 +1
     0  0 -1 +1  0
     0  0 -1  0 +1
]));

FE = SparseArrays.sparse(Int8.([
    +1 -1  0  0 +1   0  0  0  0
    +1  0 -1  0  0  +1  0  0  0
    +1  0  0 -1  0   0 +1  0  0
     0 +1 -1  0  0   0  0 +1  0
     0 +1  0 -1  0   0  0  0 +1
     0  0  0  0 +1  -1  0 +1  0
     0  0  0  0 +1   0 -1  0 +1
]));

CF = SparseArrays.sparse(Int8.([
    -1 +1  0 -1  0 +1  0
    +1  0 -1  0 +1  0 -1
     0 -1 +1 +1 -1 -1 +1
]));

model = CAGD.Model(V)
CAGD.addModelCells!(model, 1, EV)
CAGD.addModelCells!(model, 2, FE)
CAGD.addModelCells!(model, 3, CF)



atol = 1e-6;
if todisplay  displayModel(model)  end

split_model = CAGD.pairwise_decomposition(model, atol = atol)

if todisplay  displayModel(split_model)  end

congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)
#split_model = CAGD.pairwise_decomposition(congr_model, atol = atol)
#congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)

if todisplay  displayModel(congr_model)  end

gift_model, bicon_comps = CAGD.tgw(congr_model, 3, atol=atol)

if todisplay  viewExplode(gift_model)  end
