using LinearAlgebraicRepresentation
using Base.Test

function BoxCalculation(Vertices)
	Minx=minimum(Vertices[1,:])
	Maxx=maximum(Vertices[1,:])
	Miny=minimum(Vertices[2,:])
	Maxy=maximum(Vertices[2,:])
	dx=Maxx-Minx
	dy=Maxy-Miny
	Box=dx*dy
	if size(Vertices,1)==3
		Minz=minimum(Vertices[3,:])
		Maxz=maximum(Vertices[3,:])
		dz=Maxz-Minz
		Box=Box*dz
	end
	return Box
end

@testset "circle" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.circle()()[1])==4
	@test BoxCalculation(LinearAlgebraicRepresentation.circle(2., 2*pi)()[1])==16
	@test BoxCalculation(LinearAlgebraicRepresentation.circle(3, 2*pi)()[1])==36
	@test BoxCalculation(LinearAlgebraicRepresentation.circle(5, pi/2)()[1])==25
	@test size(LinearAlgebraicRepresentation.circle(3,2*pi)(60)[1],2)==60
	@test length(LinearAlgebraicRepresentation.circle(3,2*pi)(60)[2])==60
end

@testset "helix" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.helix()()[1])==8
	@test BoxCalculation(LinearAlgebraicRepresentation.helix(2, 2, 2)()[1])==64
	@test BoxCalculation(LinearAlgebraicRepresentation.helix(1, 2, 5)()[1])==40
	@test BoxCalculation(LinearAlgebraicRepresentation.helix(1, 1, 3)()[1])==12
	#=
	1, 1, 2, is contained
	in a box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LinearAlgebraicRepresentation.helix(5, 7, 9)()[1],2)==325
	@test length(LinearAlgebraicRepresentation.helix(5, 7, 9)()[2])==324
end

@testset "disk" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.disk(1., 2*pi)([36, 1])[1])==4
	@test BoxCalculation(LinearAlgebraicRepresentation.disk(2., 2*pi)()[1])==16
	@test BoxCalculation(LinearAlgebraicRepresentation.disk(2., pi)()[1])==8
	@test BoxCalculation(LinearAlgebraicRepresentation.disk(2., pi/2)()[1])==4
	@test size(LinearAlgebraicRepresentation.disk(10, pi/7)()[1],2)==75
	@test length(LinearAlgebraicRepresentation.disk(10, pi/7)()[2])==108
end

@testset "helicoid" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.helicoid()()[1])==8
	@test BoxCalculation(LinearAlgebraicRepresentation.helicoid(2,1,2,2)()[1])==64
	@test BoxCalculation(LinearAlgebraicRepresentation.helicoid(1,.5,2,5)()[1])==40
	@test BoxCalculation(LinearAlgebraicRepresentation.helicoid(1,.3,1,3)()[1])==12
	#=
	1, 1, 2, is contained in a
	box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LinearAlgebraicRepresentation.helicoid(1.3,.7,1,3)()[1],2)==327
	@test length(LinearAlgebraicRepresentation.helicoid(1.3,.7,1,3)()[2])==432
end

@testset "ring" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.ring(1.,3.,2*pi)()[1])==36
	@test BoxCalculation(LinearAlgebraicRepresentation.ring(1.,2,pi)()[1])==8
	@test BoxCalculation(LinearAlgebraicRepresentation.ring(1.,5,pi/2)()[1])==25
	#(radius*2)^2
	@test size(LinearAlgebraicRepresentation.ring(1,5,pi/2)()[1],2)==74
	@test length(LinearAlgebraicRepresentation.ring(5,10,pi/6)()[2])==36
end

@testset "cylinder" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.cylinder(1,5,2*pi)()[1])==20
	@test BoxCalculation(LinearAlgebraicRepresentation.cylinder(2,2,pi)()[1])==16
	@test BoxCalculation(LinearAlgebraicRepresentation.cylinder(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(LinearAlgebraicRepresentation.cylinder(3.4,20,pi/7)()[1],2)==74
	@test length(LinearAlgebraicRepresentation.cylinder(3.4,20,pi/7)()[2])==36
end

@testset "sphere" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.sphere(2,pi,2*pi)()[1])==64
	@test BoxCalculation(LinearAlgebraicRepresentation.sphere(6,pi,pi)()[1])==864
	@test BoxCalculation(LinearAlgebraicRepresentation.sphere(4,pi,2*pi)()[1])==8^3
	@test size(LinearAlgebraicRepresentation.sphere(2.5,pi/3,pi/5)()[1],2)==703
	@test length(LinearAlgebraicRepresentation.sphere(2.5,pi/3,pi/5)()[2])==1296
end

@testset "toroidal" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.toroidal(1,3,2*pi,2*pi)()[1])==128
	@test BoxCalculation(LinearAlgebraicRepresentation.toroidal(2,3,2*pi,2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LinearAlgebraicRepresentation.toroidal(1.3,4.6,pi/4,pi/7)()[1],2)==925
	@test length(LinearAlgebraicRepresentation.toroidal(1.3,4.6,pi/4,pi/7)()[2])==1728
end

@testset "crown" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.crown(1., 3., 2*pi)()[1])==128
	@test BoxCalculation(LinearAlgebraicRepresentation.crown(2., 3., 2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LinearAlgebraicRepresentation.crown(1.5, 5.6, pi/8)()[1],2)==481
	@test length(LinearAlgebraicRepresentation.crown(1.5, 5.6, pi/8)()[2])==864
end




@testset "ball" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.ball(1, pi, 2*pi)()[1])==8
	@test BoxCalculation(LinearAlgebraicRepresentation.ball(6, pi, pi)()[1])==864
	@test BoxCalculation(LinearAlgebraicRepresentation.ball(3, pi, 2*pi)()[1])==6^3
	@test size(LinearAlgebraicRepresentation.ball(2.6, pi/5, pi/9)()[1],2)==2813
	@test length(LinearAlgebraicRepresentation.ball(2.6, pi/5, pi/9)()[2])==2592
end

@testset "rod" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.rod(1,5,2*pi)()[1])==20
	@test BoxCalculation(LinearAlgebraicRepresentation.rod(2,2,pi)()[1])==16
	@test BoxCalculation(LinearAlgebraicRepresentation.rod(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(LinearAlgebraicRepresentation.rod(3.7,8.9,pi/9)()[1],2)==74
	@test length(LinearAlgebraicRepresentation.rod(3.7,8.9,pi/9)()[2])==1
end

@testset "hollowCyl" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowCyl(0,1.,5,2*pi)()[1])==20
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowCyl(1,2.,4,pi)()[1])==32
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowCyl(1,4,3,pi/2)()[1])==48
	#((radius*2)^2)*height
	@test size(LinearAlgebraicRepresentation.hollowCyl(3,4.,7.8,pi/5)()[1],2)==148
	@test length(LinearAlgebraicRepresentation.hollowCyl(3,4.,7.8,pi/5)()[2])==36
end

@testset "hollowBall" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowBall(1,2,pi,2*pi)([36,36,1])[1]) ==64
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowBall(1,6,pi,pi)([36,36,1])[1])== 864
	@test BoxCalculation(LinearAlgebraicRepresentation.hollowBall(2,4,pi,2*pi)([36,36,1])[1])==8^3
	#(radius*2)^3
	@test size(LinearAlgebraicRepresentation.hollowBall(1.5,6.7,pi/3,2*pi/3)()[1],2)==3700
	@test length(LinearAlgebraicRepresentation.hollowBall(1.5,6.7,pi/3,2*pi/3)()[2])==2592
end

@testset "torus" begin
	@test BoxCalculation(LinearAlgebraicRepresentation.torus(1.,2.,.5,2*pi,2*pi)()[1])==147
	@test BoxCalculation(LinearAlgebraicRepresentation.torus(2,3,.5,2*pi,2*pi)()[1])==605.0
	#(((R+r)*2)^2)*(r*2)
	@test size(LinearAlgebraicRepresentation.torus(5.2,7,.5,pi/3,pi/4)()[1],2)==4625
	@test length(LinearAlgebraicRepresentation.torus(5.2,7,.5,pi/3,pi/4)()[2])==3456
end
