
using LARLIB
using Base.Test


@testset "Integration Tests" begin

	@testset "M" begin
		@test M(0,0)==0.5
		@test M(1,0)==0.16666666666666666
		@test M(2,0)==0.08333333333333333
		@test M(3,0)==0.05
		@test M(1,1)==0.041666666666666685
		@test M(2,0)==0.08333333333333333
		@test M(2,1)==0.016666666666666663
		@test M(3,0)==0.05
		@test M(3,1)==0.008333333333333338
		@test M(3,2)==0.0023809523809523586
		@test M(3,3)==0.0008928571428571397
	end


	@testset "TT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test TT(tau, 0,0,0)==0.5
		@test TT(tau, 1,0,0)==0.16666666666666666
		@test TT(tau, 1,1,0)==0.041666666666666685
		@test TT(tau, 1,1,1)==0.0
		@test TT(tau, 2,0,0)==0.08333333333333333
		@test TT(tau, 2,1,0)==0.016666666666666663
		@test TT(tau, 2,2,0)==0.005555555555555545
		@test TT(tau, 2,2,1)==0.0
		@test TT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test II(P, 0,0,0)==0.5
		@test II(P, 1,0,0)==0.16666666666666666
		@test II(P, 0,1,0)>=0.1666666666666666
		@test II(P, 0,0,1)==0.0
		@test II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test III(P, 0,0,0)>0.166666666
		@test III(P, 0,0,0)<0.166666888
		@test III(P, 1,0,0)>0.041666666
		@test III(P, 1,0,0)<0.041666888
		@test III(P, 0,1,0)>0.041666666
		@test III(P, 0,1,0)<0.041666888
		@test III(P, 0,0,1)>0.041666666
		@test III(P, 0,0,1)<0.041666888
		@test III(P, 10,10,10)>1.3377e-11
		@test III(P, 10,10,10)<1.3388e-11
	end


	@testset "surface" begin
		V,FV = LARLIB.simplexGrid([1,1]);
		P = [V;[0 0 0 0]], FV
		@test surface(P)==1.0
		p = LARLIB.Struct([LARLIB.t(0.5,0.5,0), LARLIB.r(0,0,pi/4), P]);
		q = LARLIB.struct2lar(p);
		@test surface(q)>1.0000000
		@test surface(q)<1.0000222
	end


	@testset "volume" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test volume(P)>0.166666666
		@test volume(P)<0.166668888
	end


	@testset "firstMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test firstMoment(P)[1]<0.0416667
		@test firstMoment(P)[1]>0.0416665

		@test firstMoment(P)[2]<0.0416667
		@test firstMoment(P)[2]>0.0416665

		@test firstMoment(P)[3]<0.0416667
		@test firstMoment(P)[3]>0.0416665

		@test abs(firstMoment(P)[1]-firstMoment(P)[2])<0.00001
		@test abs(firstMoment(P)[2]-firstMoment(P)[3])<0.00001
		@test abs(firstMoment(P)[3]-firstMoment(P)[1])<0.00001
	end


	@testset "secondMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test secondMoment(P)[1]<0.0166666669
		@test secondMoment(P)[1]>0.0166666664

		@test secondMoment(P)[2]<0.0166666669
		@test secondMoment(P)[2]>0.0166666664

		@test secondMoment(P)[3]<0.0166666669
		@test secondMoment(P)[3]>0.0166666664

		@test abs(secondMoment(P)[1]-secondMoment(P)[2])<0.00001
		@test abs(secondMoment(P)[2]-secondMoment(P)[3])<0.00001
		@test abs(secondMoment(P)[3]-secondMoment(P)[1])<0.00001
	end


	@testset "inertiaProduct" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test inertiaProduct(P)[1]<0.00833666
		@test inertiaProduct(P)[1]>0.00833000

		@test inertiaProduct(P)[2]<0.00833666
		@test inertiaProduct(P)[2]>0.00833000

		@test inertiaProduct(P)[3]<0.00833666
		@test inertiaProduct(P)[3]>0.00833000
	end


	@testset "centroid" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test centroid(P)[1]<0.26
		@test centroid(P)[1]>0.24

		@test centroid(P)[2]<0.26
		@test centroid(P)[2]>0.24

		@test centroid(P)[3]<0.26
		@test centroid(P)[3]>0.24
	end


	@testset "inertiaMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test inertiaMoment(P)[1]<0.0333555
		@test inertiaMoment(P)[1]>0.0333111

		@test inertiaMoment(P)[2]<0.0333555
		@test inertiaMoment(P)[2]>0.0333111

		@test inertiaMoment(P)[3]<0.0333555
		@test inertiaMoment(P)[3]>0.0333111
	end

end

