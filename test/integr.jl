
using LARLIB
using Base.Test


@testset "Integration Tests" begin

	@testset "M" begin
		@test LARLIB.M(0,0)==0.5
		@test LARLIB.M(1,0)==0.16666666666666666
		@test LARLIB.M(2,0)==0.08333333333333333
		@test LARLIB.M(3,0)==0.05
		@test LARLIB.M(1,1)==0.041666666666666685
		@test LARLIB.M(2,0)==0.08333333333333333
		@test LARLIB.M(2,1)==0.016666666666666663
		@test LARLIB.M(3,0)==0.05
		@test LARLIB.M(3,1)==0.008333333333333338
		@test LARLIB.M(3,2)==0.0023809523809523586
		@test LARLIB.M(3,3)==0.0008928571428571397
	end


	@testset "TT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test LARLIB.TT(tau, 0,0,0)==0.5
		@test LARLIB.TT(tau, 1,0,0)==0.16666666666666666
		@test LARLIB.TT(tau, 1,1,0)==0.041666666666666685
		@test LARLIB.TT(tau, 1,1,1)==0.0
		@test LARLIB.TT(tau, 2,0,0)==0.08333333333333333
		@test LARLIB.TT(tau, 2,1,0)==0.016666666666666663
		@test LARLIB.TT(tau, 2,2,0)==0.005555555555555545
		@test LARLIB.TT(tau, 2,2,1)==0.0
		@test LARLIB.TT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test LARLIB.II(P, 0,0,0)==0.5
		@test LARLIB.II(P, 1,0,0)==0.16666666666666666
		@test LARLIB.II(P, 0,1,0)>=0.1666666666666666
		@test LARLIB.II(P, 0,0,1)==0.0
		@test LARLIB.II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.III(P, 0,0,0)>0.166666666
		@test LARLIB.III(P, 0,0,0)<0.166666888
		@test LARLIB.III(P, 1,0,0)>0.041666666
		@test LARLIB.III(P, 1,0,0)<0.041666888
		@test LARLIB.III(P, 0,1,0)>0.041666666
		@test LARLIB.III(P, 0,1,0)<0.041666888
		@test LARLIB.III(P, 0,0,1)>0.041666666
		@test LARLIB.III(P, 0,0,1)<0.041666888
		@test LARLIB.III(P, 10,10,10)>1.3377e-11
		@test LARLIB.III(P, 10,10,10)<1.3388e-11
	end


	@testset "surface" begin
		V,FV = LARLIB.simplexGrid([1,1]);
		P = [V;[0 0 0 0]], FV
		@test LARLIB.surface(P)==1.0
		p = LARLIB.Struct([LARLIB.t(0.5,0.5,0), LARLIB.r(0,0,pi/4), P]);
		q = LARLIB.struct2lar(p);
		@test LARLIB.surface(q)>1.0000000
		@test LARLIB.surface(q)<1.0000222
	end


	@testset "volume" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.volume(P)>0.166666666
		@test LARLIB.volume(P)<0.166668888
	end


	@testset "firstMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.firstMoment(P)[1]<0.0416667
		@test LARLIB.firstMoment(P)[1]>0.0416665

		@test LARLIB.firstMoment(P)[2]<0.0416667
		@test LARLIB.firstMoment(P)[2]>0.0416665

		@test LARLIB.firstMoment(P)[3]<0.0416667
		@test LARLIB.firstMoment(P)[3]>0.0416665

		@test abs(LARLIB.firstMoment(P)[1]-LARLIB.firstMoment(P)[2])<0.00001
		@test abs(LARLIB.firstMoment(P)[2]-LARLIB.firstMoment(P)[3])<0.00001
		@test abs(LARLIB.firstMoment(P)[3]-LARLIB.firstMoment(P)[1])<0.00001
	end


	@testset "secondMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.secondMoment(P)[1]<0.0166666669
		@test LARLIB.secondMoment(P)[1]>0.0166666664

		@test LARLIB.secondMoment(P)[2]<0.0166666669
		@test LARLIB.secondMoment(P)[2]>0.0166666664

		@test LARLIB.secondMoment(P)[3]<0.0166666669
		@test LARLIB.secondMoment(P)[3]>0.0166666664

		@test abs(LARLIB.secondMoment(P)[1]-LARLIB.secondMoment(P)[2])<0.00001
		@test abs(LARLIB.secondMoment(P)[2]-LARLIB.secondMoment(P)[3])<0.00001
		@test abs(LARLIB.secondMoment(P)[3]-LARLIB.secondMoment(P)[1])<0.00001
	end


	@testset "inertiaProduct" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.inertiaProduct(P)[1]<0.00833666
		@test LARLIB.inertiaProduct(P)[1]>0.00833000

		@test LARLIB.inertiaProduct(P)[2]<0.00833666
		@test LARLIB.inertiaProduct(P)[2]>0.00833000

		@test LARLIB.inertiaProduct(P)[3]<0.00833666
		@test LARLIB.inertiaProduct(P)[3]>0.00833000
	end


	@testset "centroid" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.centroid(P)[1]<0.26
		@test LARLIB.centroid(P)[1]>0.24

		@test LARLIB.centroid(P)[2]<0.26
		@test LARLIB.centroid(P)[2]>0.24

		@test LARLIB.centroid(P)[3]<0.26
		@test LARLIB.centroid(P)[3]>0.24
	end


	@testset "inertiaMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test LARLIB.inertiaMoment(P)[1]<0.0333555
		@test LARLIB.inertiaMoment(P)[1]>0.0333111

		@test LARLIB.inertiaMoment(P)[2]<0.0333555
		@test LARLIB.inertiaMoment(P)[2]>0.0333111

		@test LARLIB.inertiaMoment(P)[3]<0.0333555
		@test LARLIB.inertiaMoment(P)[3]>0.0333111
	end

end

