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
      spboundary1 = characteristicMatrix(EV)'
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
   function columninfo(col)
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
   
       
       spboundary2 = sparse(I,J,V)
       return spboundary2
   end
   

   include("./utilities.jl")
   include("./minimal_cycles.jl")
   include("./dimension_travel.jl")
   include("./planar_arrangement.jl")
   include("./spatial_arrangement.jl")
   include("./largrid.jl")
   
end
