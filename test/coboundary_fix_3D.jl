using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL, LinearAlgebra #, Plasm
GL = ViewerGL

function theshow(V,EV=[],FV=[],mysize=2.0)
    VV = [[k] for k=1:size(V,2)];
    model = (V, [VV,EV,FV])::Lar.LARmodel;

    meshes = GL.numbering(mysize)(model, GL.COLORS[1], 0.1);
    GL.VIEW(meshes);
end


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

V = [0.  2   5   7  10   2   5   3   7  3  7  0  3  3  7  0  10;
    16  16  16  16  16  13  13  11  11  8  8  5  5  2  2  0   0]
    
VV = [[v] for v in 1:size(V,2)]

size(V,2) - length(EV) + length(FV) == 2 # FV include la faccia esterna; Euler relation in d=2


#hpc0 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], VV]))
#hpc1 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], EV]))
#hpc2 = Plasm.numbering(2.5)((V,[[[k] for k=1:size(V,2)], FV]))
#Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) ])
theshow(V,EV,FV,3.0)

# patch per LAR (coboundary_1)

unsafedict = Dict()
lastitem = 0
FE = Lar.lar2cop(FV) * Lar.lar2cop(EV)' .÷ 2
edgeincidence = [sum(FE[:,k]) for k=1:size(FE,2)]
badedges = [h for h=1:length(edgeincidence) if edgeincidence[h]!=2]
lastitem = 0
for h in badedges
    edgekey = V[:,EV[h][1]], V[:,EV[h][2]]
    lastitem += 1
    unsafedict[ edgekey ] = lastitem
end
unsafedict



# building chain bases

m02 = (V,VV)
#Plasm.view(m02)
#Plasm.view([ Plasm.color("yellow")(hpc0) ])

m22 = (V,FV)::Lar.LAR
#Plasm.view(m22)
#Plasm.view([ Plasm.color("cyan")(hpc1) ])

m12 = (V,EV)::Lar.LAR
#Plasm.view(m12)
#Plasm.view([ Plasm.color("magenta")(hpc2) ])

theshow(V,EV,FV,3.0)

# building 2-skeleton

v1 = [0. 5. 10.]
c0 = [[1], [2], [3]]
c1 = [[1, 2], [2,3]]
m01 = (v1,c0)::Lar.LAR
#Plasm.view(m01)
m11 = (v1,c1)::Lar.LAR
#Plasm.view(m11)
m23v = Lar.larModelProduct( m12,m11 )
#Plasm.view(m23v)
m23o =Lar.larModelProduct( m22,m01 )
#Plasm.view(m23o)
m23 = Lar.Struct([m23o,m23v])
#Plasm.view(m23)
W,FW = Lar.struct2lar(m23)
#Plasm.view((W,FW))

theshow(W,EV,W,3.0)


# building 1-skeleton

m13o = Lar.larModelProduct( m12,m01 )
#Plasm.view(m13o)
m13v = Lar.larModelProduct( m02,m11 )
#Plasm.view(m13v)

m13 = Lar.Struct([m13o,m13v])
#Plasm.view(m13)
m13 = Lar.struct2lar(m13)
#Plasm.view(m13)

U = m13[1]
EU = m13[2]
#Plasm.view(U,EU)

UU = [[v] for v in 1:size(U,2)]
hpc0 = Plasm.numbering(2.5)((U,[[[k] for k=1:size(U,2)], UU]))
hpc1 = Plasm.numbering(2.5)((U,[[[k] for k=1:size(U,2)], EU]))
#Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) ])


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
#Plasm.view([ Plasm.color("blue")(hpc0) , Plasm.color("cyan")(hpc1) , Plasm.color("black")(hpc2)  ])


GL.VIEW( GL.numbering(3.05)((W,[WW, EW, FW]),GL.COLORS[1]) );


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
FZ[20]=[2,10,52,53]; FZ[21]=[10,18,53,54]; FZ[41]=[30,39,55,56];FZ[45]=[39,48,56,57]; 
FZ[56]=[3,11,52,53]; FZ[57]=[11,19,53,54]; FZ[58]=[28,37,55,56]; FZ[59]=[37,46,56,57]; 
#FZ[20]=[2,10,52,53]; FZ[21]=[10,18,53,54];   FZ[41]=[42,51,55,56]; FZ[45]=[39,48,56,57]; 
#FZ[56]=[3,11,52,53]; FZ[57]=[11,19,53,54];   FZ[58]=[39,48,55,56]; FZ[59]=[37,46,56,57]; 

#([2, 10, 3, 11], [10, 18, 11, 19], [37, 46, 39, 48], [39, 48, 42, 51])


ZZ = [[v] for v in 1:size(Z,2)]
zhpc0 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], ZZ]))
zhpc1 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], EZ]))
zhpc2 = Plasm.numbering(2.5)((Z,[[[k] for k=1:size(Z,2)], FZ]))
#Plasm.view([ Plasm.color("blue")(zhpc0) , Plasm.color("cyan")(zhpc1), Plasm.color("")(zhpc2) ])
#Plasm.view([ Plasm.color("blue")(zhpc0) , Plasm.color("cyan")(zhpc1), (zhpc2) ])

FZ = [[1, 2, 3, 4, 25, 26, 27, 28, 29, 30, 6, 31, 32, 33], [9, 10, 11, 12, 34, 35, 36, 37, 38, 39, 14, 40, 41, 42], [17, 18, 19, 20, 43, 44, 45, 46, 47, 48, 22, 49, 50, 51], [2, 3, 25, 26, 52], [10, 11, 34, 35, 53], [18, 19, 43, 44, 54], [4, 5, 28, 30, 6, 31, 32, 33, 7, 8], [12, 13, 37, 39, 14, 40, 41, 42, 15, 16], [20, 21, 46, 48, 22, 49, 50, 51, 23, 24], [27, 28, 29, 30, 55], [36, 37, 38, 39, 56], [45, 46, 47, 48, 57], [1, 9, 2, 10], [9, 17, 10, 18], [1, 9, 6, 14], [9, 17, 14, 22], [2, 10, 52, 53], [10, 18, 53, 54], [2, 10, 25, 34], [10, 18, 34, 43], [3, 11, 4, 12], [11, 19, 12, 20], [3, 11, 26, 35], [11, 19, 35, 44], [4, 12, 5, 13], [12, 20, 13, 21], [4, 12, 28, 37], [12, 20, 37, 46], [5, 13, 8, 16], [13, 21, 16, 24], [25, 34, 26, 35], [34, 43, 35, 44], [27, 36, 28, 37], [36, 45, 37, 46], [27, 36, 29, 38], [36, 45, 38, 47], [30, 39, 55, 56], [29, 38, 30, 39], [38, 47, 39, 48], [30, 39, 33, 42], [39, 48, 56, 57], [6, 14, 31, 40], [14, 22, 40, 49], [6, 14, 7, 15], [14, 22, 15, 23], [31, 40, 32, 41], [40, 49, 41, 50], [32, 41, 33, 42], [41, 50, 42, 51], [7, 15, 8, 16], [15, 23, 16, 24], [3, 11, 52, 53], [11, 19, 53, 54], [28, 37, 55, 56], [37, 46, 56, 57], [39,42,48,51]]

length(ZZ) - length(EZ) + length(FZ) - ((4+4) + 1) == 0  # Euler relation in d=3

GL.VIEW( GL.numbering(3.05)((Z,[ZZ, EZ,FZ]),GL.COLORS[1]) );


	cscFZ = Lar.characteristicMatrix(FZ)
	cscEZ = Lar.characteristicMatrix(EZ)
	cscFE = (cscFZ * cscEZ') .÷ 2 
    cscFE = convert(Lar.SparseArrays.SparseMatrixCSC{Int8, Int64}, cscFE);    
    Z = convert(Lar.Points, Z');

FE = Lar.cop2lar(cscFE)
for k=1:length(FE) println("$k  $(FE[k])") end

V, copEV, copFE, copCF = Lar.space_arrangement(
	Z::Lar.Points, cscEZ::Lar.ChainOp, cscFE::Lar.ChainOp);


# on error:
julia> for k=1:length(sp_idx) println("$k  $(sp_idx[k])") end
julia> for k=1:length(FZ) println("$k  $(FZ[k])") end
julia> for k=1:size(cscFE,1) println("$k  $(cscFE[k,:])") end
