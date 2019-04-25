using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
using SparseArrays

V,(VV,EV,FV,CV) = Lar.cuboid([0.5,0.5,0.5],true,[-0.5,-0.5,-0.5])
mybox = (V,CV,FV,EV)

twocubes = Lar.Struct([ mybox , Lar.t(0.3,0.4,0.5), Lar.r(pi/3,0,0), Lar.r(0,0,pi/6), mybox ])
V,CV,FV,EV = Lar.struct2lar(twocubes)
Plasm.view(V,CV)

cop_EV = Lar.coboundary_0(EV::Lar.Cells);
cop_EW = convert(Lar.ChainOp, cop_EV);
cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells);
W = convert(Lar.Points, V');

V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)


#V, copEV, copFE, multiproc = W, cop_EW, cop_FE, false
#rV, rcopEV, rcopFE = Lar.Arrangement.spatial_arrangement_1( V, copEV, copFE, multiproc )


FE = [findnz(cop_FE[k,:])[1] for k=1:size(cop_FE,1)]
FV = [collect(Set(cat(EV[e] for e in FE[f]))) for f=1:length(FE)]

i=10
Lar.Arrangement.face_int(tV, EV, FE[10, :]) 
Lar.Arrangement.face_int(tV, EV, FE[7, :])

Plasm.view(Plasm.numbering(0.25)((V,[[[k] for k=1:size(v,2)],EV,FV])))




julia> v = [-0.5 -0.5 -0.5 -0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 -0.306218 0.14641 0.10718 0.116987 0.116987 -0.383013 -0.383013 -0.110684 0.983013 0.983013 0.483013 0.483013; -0.5 -0.5 0.5 0.5 -0.5 -0.5 0.5 0.5 0.5 -0.17735 -0.263953 0.5 0.5 0.5 0.491506 -0.374519 0.924519 0.0584936 -0.17735 0.741506 -0.124519 1.17452 0.308494; -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.0910683 0.5 0.35 0.5 -0.326795 -0.326795 -0.341506 0.158494 0.408494 0.908494 0.5 0.0915064 0.591506 0.841506 1.34151]
julia> ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 9], [8, 9], [5, 7], [6, 10], [8, 10], [10, 11], [9, 11], [1, 5], [2, 6], [3, 7], [4, 12], [8, 12], [9, 13], [12, 14], [13, 14], [6, 8], [15, 16], [17, 18], [14, 15], [14, 17], [16, 19], [18, 19], [12, 19], [20, 21], [22, 23], [20, 22], [21, 23], [13, 15], [13, 20], [11, 16], [11, 21], [17, 22], [18, 23], [10, 19]]
julia> Plasm.view(Plasm.numbering(0.25)((v,[[[k] for k=1:size(v,2)],ev])))



v = [0.5 0.5 0.5 0.5 0.5 0.5 0.5; -0.5 -0.5 0.5 0.5 0.5 -0.17735 -0.263953; -0.5 0.5 -0.5 0.5 -0.0910683 0.5 0.35]
ev = Array{Int64,1}[[1, 2], [3, 5], [4, 5], [1, 3], [2, 6], [4, 6], [6, 7], [5, 7]]
sigma2 = Lar.Struct([(v,ev)])
v = [-0.5 -0.5 0.5 0.5 -0.306218; -0.5 0.5 -0.5 0.5 0.5; 0.5 0.5 0.5 0.5 0.5]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 5], [4, 5]]
sigma6 = Lar.Struct([(v,ev)])
v,ev= Lar.struct2lar(Lar.Struct([sigma2,sigma6]))
julia> Plasm.view(Plasm.numbering(0.25)((v,[[[k] for k=1:size(v,2)],ev])))



sp_idx[sigma] = [1, 2, 5, 9, 7, 8, 10, 11, 4]
sigma = 6
i = 1
v = Any[0.0 1.0 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6]]
Building batches from Hpc....
...done in 0 msec
Optimizing the octree....
   Number of input batches 75
   total number vertices    26
   Number of output batches 13
   Batch vertex media       2
...done in 0 msec
Building octree from 13 batches....
Scene number of nodes of the octree 41
Scene max depth                     4
Scene number of batches             13
...done in 0 msec
sigma = 6
i = 2
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8]]
Building batches from Hpc....
...done in 0 msec
Optimizing the octree....
   Number of input batches 105
   total number vertices    26
   Number of output batches 13
   Batch vertex media       2
...done in 0 msec
Building octree from 13 batches....
Scene number of nodes of the octree 41
Scene max depth                     4
Scene number of batches             13
...done in 0 msec
sigma = 6
i = 5
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8]]
Building batches from Hpc....
...done in 0 msec
Optimizing the octree....
   Number of input batches 105
   total number vertices    26
   Number of output batches 13
   Batch vertex media       2
...done in 0 msec
Building octree from 13 batches....
Scene number of nodes of the octree 41
Scene max depth                     4
Scene number of batches             13
...done in 0 msec
sigma = 6
i = 9
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10]]
Building batches from Hpc....
...done in 0 msec
Optimizing the octree....
   Number of input batches 131
   total number vertices    44
   Number of output batches 22
   Batch vertex media       2
...done in 0 msec
Building octree from 22 batches....
Scene number of nodes of the octree 49
Scene max depth                     4
Scene number of batches             22
...done in 0 msec
sigma = 6
i = 7
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975 0.32265 1.26603; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301 0.389316 0.116987; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10], [11, 12]]
Building batches from Hpc....
...done in 1 msec
Optimizing the octree....
   Number of input batches 164
   total number vertices    46
   Number of output batches 23
   Batch vertex media       2
...done in 0 msec
Building octree from 23 batches....
Scene number of nodes of the octree 49
Scene max depth                     4
Scene number of batches             23
...done in 0 msec
sigma = 6
i = 8
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975 0.32265 1.26603 0.533975 1.47735; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301 0.389316 0.116987 1.48301 1.21068; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14]]
Building batches from Hpc....
...done in 1 msec
Optimizing the octree....
   Number of input batches 192
   total number vertices    56
   Number of output batches 28
   Batch vertex media       2
...done in 0 msec
Building octree from 28 batches....
Scene number of nodes of the octree 57
Scene max depth                     4
Scene number of batches             28
...done in 0 msec
sigma = 6
i = 10
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975 0.32265 1.26603 0.533975 1.47735 1.26603 1.47735; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301 0.389316 0.116987 1.48301 1.21068 0.116987 0.3; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16]]
Building batches from Hpc....
...done in 1 msec
Optimizing the octree....
   Number of input batches 230
   total number vertices    62
   Number of output batches 31
   Batch vertex media       2
...done in 0 msec
Building octree from 31 batches....
Scene number of nodes of the octree 60
Scene max depth                     4
Scene number of batches             31
...done in 0 msec
sigma = 6
i = 11
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975 0.32265 1.26603 0.533975 1.47735 1.26603 1.47735 1.47735 1.47735; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301 0.389316 0.116987 1.48301 1.21068 0.116987 0.3 0.3 1.21068; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16], [17, 18]]
Building batches from Hpc....
...done in 1 msec
Optimizing the octree....
   Number of input batches 262
   total number vertices    66
   Number of output batches 33
   Batch vertex media       2
...done in 0 msec
Building octree from 33 batches....
Scene number of nodes of the octree 66
Scene max depth                     4
Scene number of batches             33
...done in 0 msec
sigma = 6
i = 4
v = Any[0.0 1.0 0.0 1.0 0.0 1.0 0.0 1.0 0.32265 0.533975 0.32265 1.26603 0.533975 1.47735 1.26603 1.47735 1.47735 1.47735 1.0 1.0; 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 1.3 1.48301 0.389316 0.116987 1.48301 1.21068 0.116987 0.3 0.3 1.21068 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
ev = Array{Int64,1}[[1, 2], [3, 4], [1, 3], [2, 4], [5, 6], [7, 8], [9, 10], [11, 12], [13, 14], [15, 16], [17, 18], [19, 20]]
Building batches from Hpc....
...done in 1 msec
Optimizing the octree....
   Number of input batches 302
   total number vertices    68
   Number of output batches 34
   Batch vertex media       2
...done in 0 msec
Building octree from 34 batches....
Scene number of nodes of the octree 67
Scene max depth                     4
Scene number of batches             34





		
v=convert(Lar.Points, V')
ev = [findnz(copEV[k,:])[1] for k=1:size(copEV,1)]
Plasm.view(Plasm.numbering(0.25)((v,[[[k] for k=1:size(v,2)],ev])))


triangulated_faces = Lar.triangulate(V, [copEV, copFE])
FVs = convert(Array{Lar.Cells}, triangulated_faces)
V = convert(Lar.Points, V')
Plasm.viewcolor(V::Lar.Points, FVs::Array{Lar.Cells})

EVs = Lar.FV2EVs(copEV, copFE) # polygonal face fragments
model = V,EVs
Plasm.view(Plasm.lar_exploded(model)(1.2,1.2,1.2))
