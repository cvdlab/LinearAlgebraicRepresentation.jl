using LinearAlgebraicRepresentation
using Test

@testset "Integration Tests" begin
	@testset "M" begin
		@test LinearAlgebraicRepresentation.M(0,0)==0.5
		@test LinearAlgebraicRepresentation.M(1,0)==0.16666666666666666
		@test LinearAlgebraicRepresentation.M(2,0)==0.08333333333333333
		@test LinearAlgebraicRepresentation.M(3,0)==0.05
		@test LinearAlgebraicRepresentation.M(1,1)==0.041666666666666685
		@test LinearAlgebraicRepresentation.M(2,0)==0.08333333333333333
		@test LinearAlgebraicRepresentation.M(2,1)==0.016666666666666663
		@test LinearAlgebraicRepresentation.M(3,0)==0.05
		@test LinearAlgebraicRepresentation.M(3,1)==0.008333333333333338
		@test LinearAlgebraicRepresentation.M(3,2)==0.0023809523809523586
		@test LinearAlgebraicRepresentation.M(3,3)==0.0008928571428571397
	end


	@testset "TT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test LinearAlgebraicRepresentation.TT(tau, 0,0,0)==0.5
		@test LinearAlgebraicRepresentation.TT(tau, 1,0,0)==0.16666666666666666
		@test LinearAlgebraicRepresentation.TT(tau, 1,1,0)==0.041666666666666685
		@test LinearAlgebraicRepresentation.TT(tau, 1,1,1)==0.0
		@test LinearAlgebraicRepresentation.TT(tau, 2,0,0)==0.08333333333333333
		@test LinearAlgebraicRepresentation.TT(tau, 2,1,0)==0.016666666666666663
		@test LinearAlgebraicRepresentation.TT(tau, 2,2,0)==0.005555555555555545
		@test LinearAlgebraicRepresentation.TT(tau, 2,2,1)==0.0
		@test LinearAlgebraicRepresentation.TT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.II(P, 0,0,0)==0.5
		@test LinearAlgebraicRepresentation.II(P, 1,0,0)==0.16666666666666666
		@test LinearAlgebraicRepresentation.II(P, 0,1,0)>=0.1666666666666666
		@test LinearAlgebraicRepresentation.II(P, 0,0,1)==0.0
		@test LinearAlgebraicRepresentation.II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.III(P, 0,0,0)>0.166666666
		@test LinearAlgebraicRepresentation.III(P, 0,0,0)<0.166666888
		@test LinearAlgebraicRepresentation.III(P, 1,0,0)>0.041666666
		@test LinearAlgebraicRepresentation.III(P, 1,0,0)<0.041666888
		@test LinearAlgebraicRepresentation.III(P, 0,1,0)>0.041666666
		@test LinearAlgebraicRepresentation.III(P, 0,1,0)<0.041666888
		@test LinearAlgebraicRepresentation.III(P, 0,0,1)>0.041666666
		@test LinearAlgebraicRepresentation.III(P, 0,0,1)<0.041666888
		@test LinearAlgebraicRepresentation.III(P, 10,10,10)>1.3377e-11
		@test LinearAlgebraicRepresentation.III(P, 10,10,10)<1.3388e-11
	end


	@testset "surface" begin
		V,FV = LinearAlgebraicRepresentation.simplexGrid([1,1]);
		P = [V;[0 0 0 0]], FV
		@test LinearAlgebraicRepresentation.surface(P)==1.0
		p = LinearAlgebraicRepresentation.Struct([LinearAlgebraicRepresentation.t(0.5,0.5,0), LinearAlgebraicRepresentation.r(0,0,pi/4), P]);
		q = LinearAlgebraicRepresentation.struct2lar(p);
		@test LinearAlgebraicRepresentation.surface(q)>1.0000000
		@test LinearAlgebraicRepresentation.surface(q)<1.0000222
	end


	@testset "volume" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.volume(P)>0.166666666
		@test LinearAlgebraicRepresentation.volume(P)<0.166668888
	end


	@testset "firstMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.firstMoment(P)[1]<0.0416667
		@test LinearAlgebraicRepresentation.firstMoment(P)[1]>0.0416665

		@test LinearAlgebraicRepresentation.firstMoment(P)[2]<0.0416667
		@test LinearAlgebraicRepresentation.firstMoment(P)[2]>0.0416665

		@test LinearAlgebraicRepresentation.firstMoment(P)[3]<0.0416667
		@test LinearAlgebraicRepresentation.firstMoment(P)[3]>0.0416665

		@test abs(LinearAlgebraicRepresentation.firstMoment(P)[1]-LinearAlgebraicRepresentation.firstMoment(P)[2])<0.00001
		@test abs(LinearAlgebraicRepresentation.firstMoment(P)[2]-LinearAlgebraicRepresentation.firstMoment(P)[3])<0.00001
		@test abs(LinearAlgebraicRepresentation.firstMoment(P)[3]-LinearAlgebraicRepresentation.firstMoment(P)[1])<0.00001
	end


	@testset "secondMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.secondMoment(P)[1]<0.0166666669
		@test LinearAlgebraicRepresentation.secondMoment(P)[1]>0.0166666664

		@test LinearAlgebraicRepresentation.secondMoment(P)[2]<0.0166666669
		@test LinearAlgebraicRepresentation.secondMoment(P)[2]>0.0166666664

		@test LinearAlgebraicRepresentation.secondMoment(P)[3]<0.0166666669
		@test LinearAlgebraicRepresentation.secondMoment(P)[3]>0.0166666664

		@test abs(LinearAlgebraicRepresentation.secondMoment(P)[1]-LinearAlgebraicRepresentation.secondMoment(P)[2])<0.00001
		@test abs(LinearAlgebraicRepresentation.secondMoment(P)[2]-LinearAlgebraicRepresentation.secondMoment(P)[3])<0.00001
		@test abs(LinearAlgebraicRepresentation.secondMoment(P)[3]-LinearAlgebraicRepresentation.secondMoment(P)[1])<0.00001
	end


	@testset "inertiaProduct" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.inertiaProduct(P)[1]<0.00833666
		@test LinearAlgebraicRepresentation.inertiaProduct(P)[1]>0.00833000

		@test LinearAlgebraicRepresentation.inertiaProduct(P)[2]<0.00833666
		@test LinearAlgebraicRepresentation.inertiaProduct(P)[2]>0.00833000

		@test LinearAlgebraicRepresentation.inertiaProduct(P)[3]<0.00833666
		@test LinearAlgebraicRepresentation.inertiaProduct(P)[3]>0.00833000
	end


	@testset "centroid" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.centroid(P)[1]<0.26
		@test LinearAlgebraicRepresentation.centroid(P)[1]>0.24

		@test LinearAlgebraicRepresentation.centroid(P)[2]<0.26
		@test LinearAlgebraicRepresentation.centroid(P)[2]>0.24

		@test LinearAlgebraicRepresentation.centroid(P)[3]<0.26
		@test LinearAlgebraicRepresentation.centroid(P)[3]>0.24
	end


	@testset "inertiaMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LinearAlgebraicRepresentation.inertiaMoment(P)[1]<0.0333555
		@test LinearAlgebraicRepresentation.inertiaMoment(P)[1]>0.0333111

		@test LinearAlgebraicRepresentation.inertiaMoment(P)[2]<0.0333555
		@test LinearAlgebraicRepresentation.inertiaMoment(P)[2]>0.0333111

		@test LinearAlgebraicRepresentation.inertiaMoment(P)[3]<0.0333555
		@test LinearAlgebraicRepresentation.inertiaMoment(P)[3]>0.0333111
	end

end

