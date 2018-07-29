
using LARLIB
using Base.Test



@testset "characteristicMatrix Tests" begin
	V,(VV,EV,FV,CV) = LARLIB.cuboid([1.,1.,1.], true); 
	@test full(LARLIB.characteristicMatrix(FV)) == [
     1  1  1  1  0  0  0  0;
	 0  0  0  0  1  1  1  1;
	 1  1  0  0  1  1  0  0;
	 0  0  1  1  0  0  1  1;
	 1  0  1  0  1  0  1  0;
	 0  1  0  1  0  1  0  1]
	@test size(LARLIB.characteristicMatrix(CV))==(1,8)
	@test typeof(LARLIB.characteristicMatrix(CV))==SparseMatrixCSC{Int8,Int64}
	@test full(LARLIB.characteristicMatrix(EV)) == [
	 1  1  0  0  0  0  0  0;
	 0  0  1  1  0  0  0  0;
	 0  0  0  0  1  1  0  0;
	 0  0  0  0  0  0  1  1;
	 1  0  1  0  0  0  0  0;
	 0  1  0  1  0  0  0  0;
	 0  0  0  0  1  0  1  0;
	 0  0  0  0  0  1  0  1;
	 1  0  0  0  1  0  0  0;
	 0  1  0  0  0  1  0  0;
	 0  0  1  0  0  0  1  0;
	 0  0  0  1  0  0  0  1]
	@test size(LARLIB.characteristicMatrix(EV))==(12,8)
	@test typeof(LARLIB.characteristicMatrix(EV))==SparseMatrixCSC{Int8,Int64}
end;
   


@testset "signed_boundary_1 Tests" begin
	V,(VV,EV,FV,CV) = LARLIB.cuboid([1.,1.,1.], true);
	signed_boundary_1 = LARLIB.boundary_1( EV::LARLIB.Cells )

	@test EV==[[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], 
		[6, 8], [1, 5], [2, 6], [3, 7], [4, 8]]
	@test length(EV)==12
	@test typeof(EV)==Array{Array{Int64,1},1}
	@test size(signed_boundary_1)==(8,12)
	@test typeof(signed_boundary_1)==SparseMatrixCSC{Int8,Int64}
	@test nnz(signed_boundary_1)==24

	@test full(LARLIB.boundary_1(EV::LARLIB.Cells))==[
	 -1   0   0   0  -1   0   0   0  -1   0   0   0;
	  1   0   0   0   0  -1   0   0   0  -1   0   0;
	  0  -1   0   0   1   0   0   0   0   0  -1   0;
	  0   1   0   0   0   1   0   0   0   0   0  -1;
	  0   0  -1   0   0   0  -1   0   1   0   0   0;
	  0   0   1   0   0   0   0  -1   0   1   0   0;
	  0   0   0  -1   0   0   1   0   0   0   1   0;
	  0   0   0   1   0   0   0   1   0   0   0   1]
end



@testset "unsigned_coboundary_1 Tests" begin
	V,(VV,EV,FV,CV) = LARLIB.cuboid([1.,1.,1.], true);
	unsigned_coboundary_1 = LARLIB.u_coboundary_1(FV,EV)

	@test size(unsigned_coboundary_1)==(6,12)
	@test nnz(unsigned_coboundary_1)==24
	@test typeof(unsigned_coboundary_1)==SparseMatrixCSC{Int8,Int64}

	@test full(unsigned_coboundary_1)==[
	 1  1  0  0  1  1  0  0  0  0  0  0;
	 0  0  1  1  0  0  1  1  0  0  0  0;
	 1  0  1  0  0  0  0  0  1  1  0  0;
	 0  1  0  1  0  0  0  0  0  0  1  1;
	 0  0  0  0  1  0  1  0  1  0  1  0;
	 0  0  0  0  0  1  0  1  0  1  0  1]
 end


  


@testset "signed_coboundary_1 Tests" begin
	signed_coboundary_1 = LARLIB.coboundary_1( FV,EV );

	@test size(signed_coboundary_1)==(6,12)
	@test typeof(signed_coboundary_1)==SparseMatrixCSC{Int8,Int64}
	@test nnz(signed_coboundary_1)==24
	@test signed_coboundary_1[1,1]==-1
	@test signed_coboundary_1[6,12]==-1
	@test full(signed_coboundary_1)==
	[-1   1   0  0   1  -1  0   0  0   0   0   0;
	  0   0  -1  1   0   0  1  -1  0   0   0   0;
	 -1   0   1  0   0   0  0   0  1  -1   0   0;
	  0  -1   0  1   0   0  0   0  0   0   1  -1;
	  0   0   0  0  -1   0  1   0  1   0  -1   0;
	  0   0   0  0   0  -1  0   1  0   1   0  -1]
end




#@testset "chaincomplex 2D Tests" begin
#	W = 
#	 [0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0  3.0  3.0  3.0  3.0;
#	  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0]
#	EW = 
#	[[1, 2],[2, 3],[3, 4],[5, 6],[6, 7],[7, 8],[9, 10],[10, 11],[11, 12],[13, 14],
#	 [14, 15],[15, 16],[1, 5],[2, 6],[3, 7],[4, 8],[5, 9],[6, 10],[7, 11],[8, 12],
#	 [9, 13],[10, 14],[11, 15],[12, 16]]
#	V,bases,coboundaries = LARLIB.chaincomplex(W,EW)
#
#	@test length(bases[1])==24	# edges
#	@test typeof(bases[1])==Array{Array{Int64,1},1}	# edges
#
#	@test length(bases[2])==9 # faces -- previously unknown !!
#	@test typeof(bases[2])==Array{Array{Int64,1},1} # faces 
#
#	@test size(coboundaries[1])==(24,16) # coboundary_1 
#	@test typeof(coboundaries[1])==SparseMatrixCSC{Int8,Int64} # coboundary_1 
#	@test nnz(coboundaries[1])==48 # coboundary_1 
#
#	@test size(coboundaries[2])==(9,24) # coboundary_2: oriented 2-cycles of faces
#	@test typeof(coboundaries[2])==SparseMatrixCSC{Int8,Int64} 
#	@test nnz(coboundaries[2]) # coboundary_2: oriented 2-cycles of faces
#	@test full(coboundaries[2]) ==
#	 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0
#	  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0
#	  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0
#	  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0
#	  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0
#	  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0
#	  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0
#	  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0  0
#	  0  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1
#end

#
#
#	
#	# Example
#	```julia
#	@test cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
#	[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
#	[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
#	
#	@test cube_2 = LARLIB.Struct([LARLIB.t(0,0,0.5), LARLIB.r(0,0,pi/3), cube_1])
#	
#	@test V,FV,EV = LARLIB.struct2lar(LARLIB.Struct([ cube_1, cube_2 ]))
#	
#	@test V,bases,coboundaries = LARLIB.chaincomplex(V,FV,EV)
#	
#	@test (EV, FV, CV), (cscEV, cscFE, cscCF) = bases,coboundaries
#
#	@test FV # bases[2]
#	18-element Array{Array{Int64,1},1}:
#	 [1, 3, 4, 6]            
#	 [2, 3, 5, 6]            
#	 [7, 8, 9, 10]           
#	 [1, 2, 3, 7, 8]         
#	 [4, 6, 9, 10, 11, 12]   
#	 [5, 6, 11, 12]          
#	 [1, 4, 7, 9]            
#	 [2, 5, 11, 13]          
#	 [2, 8, 10, 11, 13]      
#	 [2, 3, 14, 15, 16]      
#	 [11, 12, 13, 17]        
#	 [11, 12, 13, 18, 19, 20]
#	 [2, 3, 13, 17]          
#	 [2, 13, 14, 18]         
#	 [15, 16, 19, 20]        
#	 [3, 6, 12, 15, 19]      
#	 [3, 6, 12, 17]          
#	 [14, 16, 18, 20]        
#
#	@test CV # bases[3]
#	3-element Array{Array{Int64,1},1}:
#	 [2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20]
#	 [2, 3, 5, 6, 11, 12, 13, 17]                    
#	 [1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]    
#	 
#	@test cscEV # coboundaries[1]
#	34×20 SparseMatrixCSC{Int8,Int64} with 68 stored entries: ...
#
#	@test cscFE # coboundaries[2]
#	18×34 SparseMatrixCSC{Int8,Int64} with 80 stored entries: ...
#	
#	@test cscCF # coboundaries[3]
#	4×18 SparseMatrixCSC{Int8,Int64} with 36 stored entries: ...
#
#
#
#
#   
#   # Collect LAR models in a single LAR model
#   function collection2model(collection)
#      W,FW,EW = collection[1]
#      shiftV = size(W,2)
#      for k=2:length(collection)
#         V,FV,EV = collection[k]
#         W = [W V]
#         FW = [FW; FV + shiftV]
#         EW = [EW; EV + shiftV]
#         shiftV = size(W,2)
#      end
#      return W,FW,EW
#   end
#   
#   
#   
#
#	# Example
#
#	@test cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
#	[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
#	[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
#	
#	@test cube_2 = LARLIB.Struct([LARLIB.t(0,0,0.5), LARLIB.r(0,0,pi/3), cube_1])
#	
#	@test W,FW,EW = LARLIB.struct2lar(LARLIB.Struct([ cube_1, cube_2 ]))
#
#	@test V,(EV,FV,EV),(cscEV,cscFE,cscCF) = LARLIB.chaincomplex(W,FW,EW)
#
#   
#   # Map 3-cells to local bases
#	"""
#
#	"""
#   function map_3cells_to_localbases(V,CV,FV,EV,cscCF,cscFE)
#      local3cells = []
#      for c=1:length(CV)
#         cf = findnz(cscCF[c+1,:])
#         tv = triangulate(cf,V,FV,EV,cscFE,cscCF)
#         vs = sort(collect(Set(hcat(tv...))))
#         vsdict = Dict([(v,k) for (k,v) in enumerate(vs)])
#         tvs = [[vsdict[t[1]],vsdict[t[2]],vsdict[t[3]]] for t in tv]
#         v = hcat([V[:,w] for w in vs]...)
#         cell = [v,tvs]
#         append!(local3cells,[cell])
#      end
#      return local3cells
#   end
#   