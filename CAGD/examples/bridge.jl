using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
CAGD = CAGD

V = [
	0.0 2.0 4.0 1.0 2.0 3.0 1.0 3.0 0.0 4.0
	4.0 4.0 4.0 3.0 3.0 3.0 1.0 1.0 0.0 0.0
];
EV = [
	[ 1, 2], [ 2, 3], [ 1, 9], [ 2, 5],
	[ 3,10], [ 4, 5], [ 5, 6], [ 4, 7],
	[ 6, 8], [ 7, 8], [ 9,10]
];
FV = [
    [1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10],
    [            4, 5,    6, 7, 8       ],
    [1, 2,    3,                   9, 10]
]
FE = [
	[1, 2, 3, 4,  4, 5, 6, 7,  8, 9, 10, 11],
	[                   6, 7,  8, 9, 10    ],
	[1, 2, 3,        5,                  11]
]

copEV = SparseArrays.sparse(Int8.([
	1  1  0  0  0  0  0  0  0  0
	0  1  1  0  0  0  0  0  0  0
	1  0  0  0  0  0  0  0  1  0
	0  1  0  0  1  0  0  0  0  0
	0  0  1  0  0  0  0  0  0  1
	0  0  0  1  1  0  0  0  0  0
	0  0  0  0  1  1  0  0  0  0
	0  0  0  1  0  0  1  0  0  0
	0  0  0  0  0  1  0  1  0  0
	0  0  0  0  0  0  1  1  0  0
	0  0  0  0  0  0  0  0  1  1
]))

copFE = SparseArrays.sparse(Int8.([
	1  1  1  2  1  1  1  1  1  1  1
	0  0  0  0  0  1  1  1  1  1  0
	1  1  1  0  1  0  0  0  0  0  1
]))

copFV = SparseArrays.sparse(Int8.([
    1  2  1  1  2  1  1  1  1  1
    0  0  0  1  1  1  1  1  0  0
    1  1  1  0  0  0  0  0  1  1
]))

#--- Using CAGD
model = CAGD.Model(V)
CAGD.addModelCells!(model, 1, EV, signed = false)
CAGD.addModelCells!(model, 2, FE, signed = false)
copEV = model.T[1]
copFE = model.T[2]
#--- Otherwise



copFV = map(x -> floor(Int8, x / 2), copFE * copEV)

copFE = map(x -> floor(Int8, x / 2), copFV * copEV')

# EXAMPLE 2
#

V = [
    1.0 3.0  0.0 1.0 3.0 4.0  2.0  2.0  1.0 3.0
    4.0 4.0  3.0 3.0 3.0 3.0  2.0  1.0  0.0 0.0
]

copEV = SparseArrays.sparse(Int8.([
    1 0 1 0 0  0 0 0 0 0
    1 0 0 1 0  0 0 0 0 0
    0 1 0 0 1  0 0 0 0 0
    0 1 0 0 0  1 0 0 0 0
    0 0 1 1 0  0 0 0 0 0
    0 0 0 0 1  1 0 0 0 0
    0 0 0 1 0  0 1 0 0 0
    0 0 0 0 1  0 1 0 0 0
    0 0 0 0 0  0 1 1 0 0
    0 0 0 0 0  0 0 1 1 0
    0 0 0 0 0  0 0 1 0 1
    0 0 0 0 0  0 0 0 1 1
]))

copFE = SparseArrays.sparse(Int8.([
    1 1 0 0 1  0 0 0 0 0  0 0
    0 0 1 1 0  1 0 0 0 0  0 0
    0 0 0 0 0  0 0 0 0 1  1 1
    1 1 1 1 1  1 2 2 2 1  1 1
]))
