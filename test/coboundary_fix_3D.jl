using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL, LinearAlgebra, Plasm
GL = ViewerGL

FV = [[1,2,3,4,5,17,16,12],
[1,2,3,4,6,7,8,9,10,11,12,13,14,15],
[4,5,9,11,12,13,14,15,16,17],
[2,3,6,7], [8,9,10,11]]

FV = sort(map(sort,FV)) # forma canonica

FE = [[1,2,3,4,9,20,17,5],
[1,6,10,7,3,8,11,12,14,15,19,18,16,5],
[4,9,20,17,16,18,19,15,13,8],
[2,10,6,7], [11,12,13,14]]

FE = sort(map(sort,FE)) # forma canonica

EV = [[1,2],[2,3],[3,4],[4,5],[1,12],[2,6],[3,7],[4,9],[5,17],[6,7],[8,9],
[8,10],[9,11],[10,11],[11,15],[12,13],[12,16],[13,14],[14,15],[16,17]]

EV = sort(map(sort,EV)) # forma canonica

V = [0   2   5   7  10   2   5   3   7  3  7  0  3  3  7  0  10;
    16  16  16  16  16  13  13  11  11  8  8  5  5  2  2  0   0]
    
VV = [[v] for v in 1:size(V,2)]

hpc0 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], VV]))
hpc1 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], EV]))
hpc2 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], FV]))
Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) ])


# building chain bases

m02 = (V,VV)
Plasm.view(m02)
Plasm.view([ Plasm.color("yellow")(hpc0) ])

m22 = (V,FV)::Lar.LAR
Plasm.view(m22)
Plasm.view([ Plasm.color("cyan")(hpc1) ])

m12 = (V,EV)::Lar.LAR
Plasm.view(m12)
Plasm.view([ Plasm.color("magenta")(hpc2) ])


# building 2-skeleton

v1 = [0. 5. 10.]
c0 = [[1], [2], [3]]
c1 = [[1, 2], [2,3]]
m01 = (v1,c0)::Lar.LAR
Plasm.view(m01)

m11 = (v1,c1)::Lar.LAR
Plasm.view(m11)

m23v = Lar.larModelProduct( m12,m11 )
Plasm.view(m23v)
m23o =Lar.larModelProduct( m22,m01 )
Plasm.view(m23o)

m23 = Lar.Struct([m23o,m23v])
Plasm.view(m23)
W,FW = Lar.struct2lar(m23)
Plasm.view((W,FW))


# building 1-skeleton

m13o = Lar.larModelProduct( m12,m01 )
Plasm.view(m13o)
m13v = Lar.larModelProduct( m02,m11 )
Plasm.view(m13v)

m13 = Lar.Struct([m13o,m13v])
Plasm.view(m13)
m13 = Lar.struct2lar(m13)
Plasm.view(m13)

U = m13[1]
EU = m13[2]
Plasm.view(U,EU)

UU = [[v] for v in 1:size(U,2)]
hpc0 = Plasm.numbering(2.5)((U,[[[k] for k=1:size(U,2)], UU]))
hpc1 = Plasm.numbering(2.5)((U,[[[k] for k=1:size(U,2)], EU]))
Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) ])


wdict = Dict{Vector{Float64},Int64}()
for k=1:size(W,2) 
    wdict[W[:,k]] = k 
end
EW = copy(EU)
for h =1:length(EU) 
    for k=1:length(EU[h])
        oldindex = EU[h][k]
        newindex = wdict[ U[:,oldindex] ]; 
        EW[h][k] = newindex
    end 
end
EW = sort(map(sort, EW))

WW = [[v] for v in 1:size(U,2)]
hpc0 = Plasm.numbering(2.5)((W,[[[k] for k=1:size(W,2)], WW]))
hpc1 = Plasm.numbering(2.5)((W,[[[k] for k=1:size(W,2)], EW]))
hpc2 = Plasm.numbering(2.5)((W,[[[k] for k=1:size(W,2)], FW]))
Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) , Plasm.color("black")(hpc2)  ])


julia> FW[20],FW[21],  FW[41],FW[45]
([2, 10, 3, 11], [10, 18, 11, 19], [37, 46, 39, 48], [39, 48, 42, 51])

julia> EW[4],EW[24],EW[43],  EW[59],EW[76],EW[90]
([2, 3], [10, 11], [18, 19], [28, 30], [37, 39], [46, 48])

W
# 3×51 Matrix{Float64}:
Z = [W zeros(3,6)]
# 3×57 Matrix{Float64}:
EZ = [EW; [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]]
# 100-element Vector{Vector{Int64}}:


Z[:,52] = (W[:,2]+W[:,3])./ 2; Z[:,53] = (W[:,10]+W[:,11])./ 2; Z[:,54] = (W[:,18]+W[:,19])./ 2; 
Z[:,55] = (W[:,28]+W[:,30])./ 2; Z[:,56] = (W[:,37]+W[:,39])./ 2; Z[:,57] = (W[:,46]+W[:,48])./ 2; 

EZ[4]=[2,52]; EZ[24]=[10,53]; EZ[43]=[18,54]; EZ[59]=[28,55]; EZ[76]=[37,56]; EZ[90]=[46,57]; 
EZ[95]=[3,52]; EZ[96]=[11,53]; EZ[97]=[19,54]; EZ[98]=[30,55]; EZ[99]=[39,56]; EZ[100]=[48,57]; 
EZ = [EZ; [[0,0],[0,0],[0,0],[0,0]]]
EZ[101]=[52,53]; EZ[102]=[53,54];      EZ[103]=[55,56]; EZ[104]=[56,57];  

FZ = [FW;  [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]]
FZ[20]=[2,10,52,53]; FZ[21]=[10,18,53,54];   FZ[41]=[30,39,55,56]; FZ[45]=[39,48,56,57]; 
FZ[56]=[3,11,52,53]; FZ[57]=[11,19,53,54];   FZ[58]=[28,37,55,56]; FZ[59]=[37,46,56,57]; 

#([2, 10, 3, 11], [10, 18, 11, 19], [37, 46, 39, 48], [39, 48, 42, 51])


ZZ = [[v] for v in 1:size(Z,2)]
zhpc0 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], ZZ]))
zhpc1 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], EZ]))
zhpc2 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], FZ]))
#Plasm.view([ Plasm.color("blue")(zhpc0) , Plasm.color("cyan")(zhpc1), Plasm.color("")(zhpc2) ])
Plasm.view([ Plasm.color("blue")(zhpc0) , Plasm.color("cyan")(zhpc1), (zhpc2) ])



	cscFZ = Lar.characteristicMatrix(FZ)
	cscEZ = Lar.characteristicMatrix(EZ)
	cscFE = (cscFZ * cscEZ') .÷ 2 
    cscFE = convert(Lar.SparseArrays.SparseMatrixCSC{Int8, Int64}, cscFE);    
    Z = convert(Lar.Points, Z');

larFE = Lar.cop2lar(cscFE)
for k=1:length(larFE) println("$k  $(larFE[k])") end

V, copEV, copFE, copCF = Lar.space_arrangement(
	Z::Lar.Points, cscEZ::Lar.ChainOp, cscFE::Lar.ChainOp);


# on error:
julia> for k=1:length(sp_idx) println("$k  $(sp_idx[k])") end
julia> for k=1:length(FZ) println("$k  $(FZ[k])") end
