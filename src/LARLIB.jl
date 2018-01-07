module LARLIB

   using NearestNeighbors
   using DataStructures
   using NearestNeighbors
   using IntervalTrees
   using TRIANGLE
   
   const Verts = Array{Float64, 2}
   const Cells = SparseMatrixCSC{Int8, Int}
   const Cell = SparseVector{Int8, Int}
   const LarCells = Array{Array{Int, 1}, 1}
   
   
   # Characteristic matrix $M_2$, i.e. M(FV)
   function characteristicMatrix(FV)
      I,J,V = Int64[],Int64[],Int8[] 
      for f=1:length(FV)
         for k in FV[f]
            push!(I,f)
            push!(J,k)
            push!(V,1)
         end
      end
      M_2 = sparse(I,J,V)
      return M_2
   end
   
   
   # Computation of sparse boundary $C_1 \to C_0$
   function boundary1(EV)
      spboundary1 = LARLIB.characteristicMatrix(EV)'
      for e = 1:length(EV)
         spboundary1[EV[e][1],e] = -1
      end
      return spboundary1
   end
   
   
   # Computation of sparse uboundary2
   function uboundary2(FV,EV)
      cscFV = characteristicMatrix(FV)
      cscEV = characteristicMatrix(EV)
      temp = cscFV * cscEV'
      I,J,V = Int64[],Int64[],Int8[]
      for j=1:size(temp,2)
         for i=1:size(temp,1)
            if temp[i,j] == 2
               push!(I,i)
               push!(J,j)
               push!(V,1)
            end
         end
      end
      sp_uboundary2 = sparse(I,J,V)
      return sp_uboundary2
   end
   
   
   
   # Local storage
   function columninfo(infos,EV,next,col)
       infos[1,col] = 1
       infos[2,col] = next
       infos[3,col] = EV[next][1]
       infos[4,col] = EV[next][2]
       vpivot = infos[4,col]
   end
   
   
   # Initialization
   function boundary2(FV,EV)
       sp_u_boundary2 = uboundary2(FV,EV)
       larEV = characteristicMatrix(EV)
       # unsigned incidence relation
       FE = [findn(sp_u_boundary2[f,:]) for f=1:size(sp_u_boundary2,1) ]
       I,J,V = Int64[],Int64[],Int8[]
       vedges = [findn(larEV[:,v]) for v=1:size(larEV,2)]
   
       # Loop on faces
       for f=1:length(FE)
           fedges = Set(FE[f])
           next = pop!(fedges)
           col = 1
           infos = zeros(Int64,(4,length(FE[f])))
           vpivot = infos[4,col]
           vpivot = columninfo(infos,EV,next,col)
           while fedges != Set()
               nextedge = intersect(fedges, Set(vedges[vpivot]))
               fedges = setdiff(fedges,nextedge)
               next = pop!(nextedge)
               col += 1
               vpivot = columninfo(infos,EV,next,col)
               if vpivot == infos[4,col-1]
                   infos[3,col],infos[4,col] = infos[4,col],infos[3,col]
                   infos[1,col] = -1
                   vpivot = infos[4,col]
               end
           end
           for j=1:size(infos,2)
               push!(I, f)
               push!(J, infos[2,j])
               push!(V, infos[1,j])
           end
       end
       
       spboundary2 = sparse(I,J,V)
       return spboundary2
   end
   
   # Chain 3-complex construction
   function chaincomplex(W,FW,EW)
       V = convert(Array{Float64,2},W')
       EV = characteristicMatrix(EW)
       FE = boundary2(FW,EW)
       V,cscEV,cscFE,cscCF = LARLIB.spatial_arrangement(V,EV,FE)
       ne,nv = size(cscEV)
       nf = size(cscFE,1)
       nc = size(cscCF,1)
       EV = [findn(cscEV[e,:]) for e=1:ne]
       FV = [collect(Set(vcat([EV[e] for e in findn(cscFE[f,:])]...)))  for f=1:nf]
       CV = [collect(Set(vcat([FV[f] for f in findn(cscCF[c,:])]...)))  for c=2:nc]
       function ord(cells)
           return [sort(cell) for cell in cells]
       end
       temp = copy(cscEV')
       for k=1:size(temp,2)
           h = findn(temp[:,k])[1]
           temp[h,k] = -1
       end    
       cscEV = temp'
       bases, coboundaries = (ord(EV),ord(FV),ord(CV)), (cscEV,cscFE,cscCF)
       return V',bases,coboundaries
   end
   

   include("./utilities.jl")
   include("./minimal_cycles.jl")
   include("./dimension_travel.jl")
   include("./planar_arrangement.jl")
   include("./spatial_arrangement.jl")
   include("./largrid.jl")
   
end
