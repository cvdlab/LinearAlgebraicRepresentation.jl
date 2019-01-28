using Test
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

@testset "Integration Tests" begin
	@testset "M" begin
		@test Lar.M(0,0)==0.5
		@test Lar.M(1,0)==0.16666666666666666
		@test Lar.M(2,0)==0.08333333333333333
		@test Lar.M(3,0)==0.05
		@test Lar.M(1,1)==0.041666666666666685
		@test Lar.M(2,0)==0.08333333333333333
		@test Lar.M(2,1)==0.016666666666666663
		@test Lar.M(3,0)==0.05
		@test Lar.M(3,1)==0.008333333333333338
		@test Lar.M(3,2)==0.0023809523809523586
		@test Lar.M(3,3)==0.0008928571428571397
	end


	@testset "TT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test Lar.TT(tau, 0,0,0)==0.5
		@test Lar.TT(tau, 1,0,0)==0.16666666666666666
		@test Lar.TT(tau, 1,1,0)==0.041666666666666685
		@test Lar.TT(tau, 1,1,1)==0.0
		@test Lar.TT(tau, 2,0,0)==0.08333333333333333
		@test Lar.TT(tau, 2,1,0)==0.016666666666666663
		@test Lar.TT(tau, 2,2,0)==0.005555555555555545
		@test Lar.TT(tau, 2,2,1)==0.0
		@test Lar.TT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test Lar.II(P, 0,0,0)==0.5
		@test Lar.II(P, 1,0,0)==0.16666666666666666
		@test Lar.II(P, 0,1,0)>=0.1666666666666666
		@test Lar.II(P, 0,0,1)==0.0
		@test Lar.II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.III(P, 0,0,0)>0.166666666
		@test Lar.III(P, 0,0,0)<0.166666888
		@test Lar.III(P, 1,0,0)>0.041666666
		@test Lar.III(P, 1,0,0)<0.041666888
		@test Lar.III(P, 0,1,0)>0.041666666
		@test Lar.III(P, 0,1,0)<0.041666888
		@test Lar.III(P, 0,0,1)>0.041666666
		@test Lar.III(P, 0,0,1)<0.041666888
		@test Lar.III(P, 10,10,10)>1.3377e-11
		@test Lar.III(P, 10,10,10)<1.3388e-11
	end


	@testset "surface" begin
		V,FV = Lar.simplexGrid([1,1]);
		P = [V;[0 0 0 0]], FV
		@test Lar.surface(P)==1.0
		p = Lar.Struct([Lar.t(0.5,0.5,0), Lar.r(0,0,pi/4), P]);
		q = Lar.struct2lar(p);
		@test Lar.surface(q)>1.0000000
		@test Lar.surface(q)<1.0000222
	end


	@testset "volume" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.volume(P)>0.166666666
		@test Lar.volume(P)<0.166668888
	end


	@testset "firstMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.firstMoment(P)[1]<0.0416667
		@test Lar.firstMoment(P)[1]>0.0416665

		@test Lar.firstMoment(P)[2]<0.0416667
		@test Lar.firstMoment(P)[2]>0.0416665

		@test Lar.firstMoment(P)[3]<0.0416667
		@test Lar.firstMoment(P)[3]>0.0416665

		@test abs(Lar.firstMoment(P)[1]-Lar.firstMoment(P)[2])<0.00001
		@test abs(Lar.firstMoment(P)[2]-Lar.firstMoment(P)[3])<0.00001
		@test abs(Lar.firstMoment(P)[3]-Lar.firstMoment(P)[1])<0.00001
	end


	@testset "secondMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.secondMoment(P)[1]<0.0166666669
		@test Lar.secondMoment(P)[1]>0.0166666664

		@test Lar.secondMoment(P)[2]<0.0166666669
		@test Lar.secondMoment(P)[2]>0.0166666664

		@test Lar.secondMoment(P)[3]<0.0166666669
		@test Lar.secondMoment(P)[3]>0.0166666664

		@test abs(Lar.secondMoment(P)[1]-Lar.secondMoment(P)[2])<0.00001
		@test abs(Lar.secondMoment(P)[2]-Lar.secondMoment(P)[3])<0.00001
		@test abs(Lar.secondMoment(P)[3]-Lar.secondMoment(P)[1])<0.00001
	end


	@testset "inertiaProduct" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.inertiaProduct(P)[1]<0.00833666
		@test Lar.inertiaProduct(P)[1]>0.00833000

		@test Lar.inertiaProduct(P)[2]<0.00833666
		@test Lar.inertiaProduct(P)[2]>0.00833000

		@test Lar.inertiaProduct(P)[3]<0.00833666
		@test Lar.inertiaProduct(P)[3]>0.00833000
	end


	@testset "centroid" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.centroid(P)[1]<0.26
		@test Lar.centroid(P)[1]>0.24

		@test Lar.centroid(P)[2]<0.26
		@test Lar.centroid(P)[2]>0.24

		@test Lar.centroid(P)[3]<0.26
		@test Lar.centroid(P)[3]>0.24
	end


	@testset "inertiaMoment" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test Lar.inertiaMoment(P)[1]<0.0333555
		@test Lar.inertiaMoment(P)[1]>0.0333111

		@test Lar.inertiaMoment(P)[2]<0.0333555
		@test Lar.inertiaMoment(P)[2]>0.0333111

		@test Lar.inertiaMoment(P)[3]<0.0333555
		@test Lar.inertiaMoment(P)[3]>0.0333111
	end

end

