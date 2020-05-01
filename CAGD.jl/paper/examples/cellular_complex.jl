# Input data for generation of the cellular complex shown in figure of Data 2.1.2

V = [
    0.0 1.5 3.0 1.0 1.5 2.0 1.0 1.5 2.0 0.0 1.5 3.0;
    0.0 0.0 0.0 1.0 1.0 1.0 2.0 2.0 2.0 3.0 3.0 3.0
];
EV = [
    [1,  2], [2,  3], [ 4,  5], [ 5,  6],
    [7,  8], [8,  9], [10, 11], [11, 12],
    [1, 10], [4,  7], [ 6,  9], [ 3, 12],
    [2,  5], [8, 11]
];
FV = [
    [1, 2, 4, 5, 7, 8, 10, 11],
    [2, 3, 5, 6, 8, 9, 11, 12],
    [4, 5, 6, 7, 8, 9]
];

# The function K returns the Julia SparseArrays matrix providing
#  the characteristic matrix K𝑟 given as input an array CV
#  specifying each 𝑟-cell as array of vertex indices.

function K( CV )
    I = vcat( [ [k for h in CV[k]] for k=1:length(CV) ]...)
    J = vcat( CV...)
    X = Int8[1 for k=1:length(I)]
    return SparseArrays.sparse(I, J, X)
end

# Generation of the sparse binary matrices

copEV = K(EV);
copFV = K(FV);

Matrix(copEV)
Matrix(copFV)

if todisplay
    VV = [[k] for k in 1:size(V,2)];
    GL.VIEW( GL.numbering(.5)( (V, [VV, EV]),GL.COLORS[1],0.1 ) );
end

# Boundary matrix generation
# we construct the sparse matrix [𝜕2] from the product of characteristic matrices EV and the transposed FV

Matrix( copEV * copFV' )

copEF = (copEV * copFV') .÷ 2;

Matrix(copEF)

b1 = (copEF * sparse([1 1 1])') .% 2
b2 = (copEF * sparse([1 1 0])') .% 2

f1 = SparseArrays.findnz(b1)[1]
f2 = SparseArrays.findnz(b2)[1]

e1 = EV[f1]
e2 = EV[f2]

if todisplay
    GL.VIEW( GL.numbering(.5)( (V, [VV, e1]),GL.COLORS[1],0.1 ) );
    GL.VIEW( GL.numbering(.5)( (V, [VV, e2]),GL.COLORS[1],0.1 ) );
end