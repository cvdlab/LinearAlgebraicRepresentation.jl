using Test
using SparseArrays
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

@testset "interface.jl file Tests" begin
	@testset "characteristicMatrix Tests" begin
		V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true); 
		@test Matrix(Lar.characteristicMatrix(FV)) == [
		 1  1  1  1  0  0  0  0;
		 0  0  0  0  1  1  1  1;
		 1  1  0  0  1  1  0  0;
		 0  0  1  1  0  0  1  1;
		 1  0  1  0  1  0  1  0;
		 0  1  0  1  0  1  0  1]
		@test size(Lar.characteristicMatrix(CV))==(1,8)
		@test typeof(Lar.characteristicMatrix(CV))==SparseMatrixCSC{Int8,Int64}
		@test Matrix(Lar.characteristicMatrix(EV)) == [
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
		@test size(Lar.characteristicMatrix(EV))==(12,8)
		@test typeof(Lar.characteristicMatrix(EV))==SparseMatrixCSC{Int8,Int64}
	end

	@testset "signed_boundary_1 Tests" begin
		V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true);
		signed_boundary_1 = Lar.boundary_1( EV::Lar.Cells )

		@test EV==[[1, 2], [3, 4], [5, 6], [7, 8], [1, 3], [2, 4], [5, 7], 
			[6, 8], [1, 5], [2, 6], [3, 7], [4, 8]]
		@test length(EV)==12
		@test typeof(EV)==Array{Array{Int64,1},1}
		@test size(signed_boundary_1)==(8,12)
		@test typeof(signed_boundary_1)==SparseMatrixCSC{Int8,Int64}
		@test nnz(signed_boundary_1)==24

		@test Matrix(Lar.boundary_1(EV::Lar.Cells))==[
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
		V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true);
		unsigned_coboundary_1 = Lar.u_coboundary_1(FV,EV)

		@test size(unsigned_coboundary_1)==(6,12)
		@test nnz(unsigned_coboundary_1)==24
		@test typeof(unsigned_coboundary_1)==SparseMatrixCSC{Int8,Int64}

		@test Matrix(unsigned_coboundary_1)==[
		 1  1  0  0  1  1  0  0  0  0  0  0;
		 0  0  1  1  0  0  1  1  0  0  0  0;
		 1  0  1  0  0  0  0  0  1  1  0  0;
		 0  1  0  1  0  0  0  0  0  0  1  1;
		 0  0  0  0  1  0  1  0  1  0  1  0;
		 0  0  0  0  0  1  0  1  0  1  0  1]
	 end

	@testset "signed_coboundary_1 Tests" begin
		V,(VV,EV,FV,CV) = Lar.cuboid([1.,1.,1.], true);
		signed_coboundary_1 = Lar.coboundary_1( V,FV,EV );

		@test size(signed_coboundary_1)==(6,12)
		@test typeof(signed_coboundary_1)==SparseArrays.SparseMatrixCSC{Int8,Int64}
		@test nnz(signed_coboundary_1)==24
		@test signed_coboundary_1[1,1]==1
		@test signed_coboundary_1[6,12]==1
		@test Matrix(signed_coboundary_1)==
			[1  -1   0   0  -1  1   0   0   0   0   0  0;
			 0   0   1  -1   0  0  -1   1   0   0   0  0;
			 1   0  -1   0   0  0   0   0  -1   1   0  0;
			 0   1   0  -1   0  0   0   0   0   0  -1  1;
			 0   0   0   0   1  0  -1   0  -1   0   1  0;
			 0   0   0   0   0  1   0  -1   0  -1   0  1]
	end

@testset "chaincomplex 2D Tests" begin
	W = 
	 [0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  2.0  2.0  2.0  2.0  3.0  3.0  3.0  3.0;
	  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0  0.0  1.0  2.0  3.0]
	EW = 
	[[1, 2],[2, 3],[3, 4],[5, 6],[6, 7],[7, 8],[9, 10],[10, 11],[11, 12],[13, 14],
	 [14, 15],[15, 16],[1, 5],[2, 6],[3, 7],[4, 8],[5, 9],[6, 10],[7, 11],[8, 12],
	 [9, 13],[10, 14],[11, 15],[12, 16]]
	V,bases,coboundaries = Lar.chaincomplex(W,EW)

	@test length(bases[1])==24	# edges
	@test typeof(bases[1])==Array{Array{Int64,1},1}	# edges

	@test length(bases[2])==9 # faces -- previously unknown !!
	@test typeof(bases[2])==Array{Array{Int64,1},1} # faces 

	@test size(coboundaries[1])==(24, 16) # coboundary_1 
	@test typeof(coboundaries[1])==SparseArrays.SparseMatrixCSC{Int8,Int64} 
	# coboundary_1 
	@test SparseArrays.nnz(coboundaries[1])==48 # coboundary_1 

	@test size(coboundaries[2])==(9,24) # coboundary_2: oriented 2-cycles of faces
	@test typeof(coboundaries[2])==SparseArrays.SparseMatrixCSC{Int8,Int64} 
	@test SparseArrays.nnz(coboundaries[2])==36 
	# coboundary_2: oriented 2-cycles of faces
	@test Matrix(coboundaries[2]) ==
	[-1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0;
	  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0;
	  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0;
	  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0;
	  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0;
	  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0;
	  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0;
	  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1  0  0;
	  0  0  0  0  0  0  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  1 -1]
end

#	@testset "chaincomplex 3D Tests" begin
#		cube_1 = ([0 0 0 0 1 1 1 1; 0 0 1 1 0 0 1 1; 0 1 0 1 0 1 0 1], 
#		[[1,2,3,4],[5,6,7,8],[1,2,5,6],[3,4,7,8],[1,3,5,7],[2,4,6,8]], 
#		[[1,2],[3,4],[5,6],[7,8],[1,3],[2,4],[5,7],[6,8],[1,5],[2,6],[3,7],[4,8]] )
#	
#		cube_2 = Lar.Struct(
#			[Lar.t(
#			0,0,0.5), Lar.r(0,0,pi/3), cube_1])
#		V,FV,EV = Lar.struct2lar(
#			Lar.Struct([ cube_1, cube_2 ]))
#			
#		V,bases,coboundaries = Lar.chaincomplex(V,FV,EV)
#		(EV, FV, CV), (cscEV, cscFE, cscCF) = bases,coboundaries
#
#		@test typeof(FV)==Array{Array{Int64,1},1} # bases[2]
#		@test length(FV)==18 # bases[2]
#		@test FV==[[1,3,4,6],[2,3,5,6],[7,8,9,10],[1,2,3,7,8],[4,6,9,10,11,12],[5,6,11,12],
#		[1,4,7,9],[2,5,11,13],[2,8,10,11,13],[2,3,14,15,16],[11,12,13,17],
#		[11,12,13,18,19,20],[2,3,13,17],[2,13,14,18],[15,16,19,20],[3,6,12,15,19],
#		[3,6,12,17],[14,16,18,20]]
#		@test typeof(CV)==Array{Array{Int64,1},1} # bases[3]
#		@test length(CV)==3 # bases[3]
#		@test CV==[[2, 3, 5, 6, 11, 12, 13, 14, 15, 16, 18, 19, 20],
#		[2, 3, 5, 6, 11, 12, 13, 17],                
#		[1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17]]   
#	 
#		@test typeof(cscEV)==SparseArrays.SparseMatrixCSC{Int8,Int64} # coboundaries[1]
#		@test size(cscEV)==(34,20) # coboundaries[1]
#		@test nnz(cscEV)==68 # coboundaries[1]
#
#		@test typeof(cscFE)==SparseArrays.SparseMatrixCSC{Int8,Int64} # coboundaries[2]
#		@test size(cscFE)==(18,34) # coboundaries[2]
#		@test nnz(cscFE)==80 # coboundaries[2]
#	
#		@test typeof(cscCF)==SparseArrays.SparseMatrixCSC{Int8,Int64} # coboundaries[3]
#		@test size(cscCF)==(4,18) # coboundaries[3]
#		@test nnz(cscCF)==36 # coboundaries[3]
#	end
end