using LARLIB

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
	@test BoxCalculation(LARLIB.circle()()[1])==4
	@test BoxCalculation(LARLIB.circle(2., 2*pi)()[1])==16
	@test BoxCalculation(LARLIB.circle(3, 2*pi)()[1])==36
	@test BoxCalculation(LARLIB.circle(5, pi/2)()[1])==25
	@test size(LARLIB.circle(3,2*pi)(60)[1],2)==60
	@test length(LARLIB.circle(3,2*pi)(60)[2])==60
end



@testset "helix" begin
	@test BoxCalculation(LARLIB.helix()()[1])==8
	@test BoxCalculation(LARLIB.helix(2, 2, 2)()[1])==64
	@test BoxCalculation(LARLIB.helix(1, 2, 5)()[1])==40
	@test BoxCalculation(LARLIB.helix(1, 1, 3)()[1])==12
	#=
	1, 1, 2, is contained
	in a box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LARLIB.helix(5, 7, 9)()[1],2)==325
	@test length(LARLIB.helix(5, 7, 9)()[2])==324
end



@testset "disk" begin
	@test BoxCalculation(LARLIB.disk(1., 2*pi)([36, 1])[1])==4
	@test BoxCalculation(LARLIB.disk(2., 2*pi)()[1])==16
	@test BoxCalculation(LARLIB.disk(2., pi)()[1])==8
	@test BoxCalculation(LARLIB.disk(2., pi/2)()[1])==4
	@test size(LARLIB.disk(10, pi/7)()[1],2)==75
	@test length(LARLIB.disk(10, pi/7)()[2])==108
end



@testset "helicoid" begin
	@test BoxCalculation(LARLIB.helicoid()()[1])==8
	@test BoxCalculation(LARLIB.helicoid(2,1,2,2)()[1])==64
	@test BoxCalculation(LARLIB.helicoid(1,.5,2,5)()[1])==40
	@test BoxCalculation(LARLIB.helicoid(1,.3,1,3)()[1])==12
	#=
	1, 1, 2, is contained in a
	box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LARLIB.helicoid(1.3,.7,1,3)()[1],2)==327
	@test length(LARLIB.helicoid(1.3,.7,1,3)()[2])==432
end


@testset "ring" begin
	@test BoxCalculation(LARLIB.ring(1.,3.,2*pi)()[1])==36
	@test BoxCalculation(LARLIB.ring(1.,2,pi)()[1])==8
	@test BoxCalculation(LARLIB.ring(1.,5,pi/2)()[1])==25
	#(radius*2)^2
	@test size(LARLIB.ring(1,5,pi/2)()[1],2)==74
	@test length(LARLIB.ring(5,10,pi/6)()[2])==36
end



@testset "cylinder" begin
	@test BoxCalculation(LARLIB.cylinder(1,5,2*pi)()[1])==20
	@test BoxCalculation(LARLIB.cylinder(2,2,pi)()[1])==16
	@test BoxCalculation(LARLIB.cylinder(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(LARLIB.cylinder(3.4,20,pi/7)()[1],2)==74
	@test length(LARLIB.cylinder(3.4,20,pi/7)()[2])==36
end



@testset "sphere" begin
	@test BoxCalculation(LARLIB.sphere(2,pi,2*pi)()[1])==64
	@test BoxCalculation(LARLIB.sphere(6,pi,pi)()[1])==864
	@test BoxCalculation(LARLIB.sphere(4,pi,2*pi)()[1])==8^3
	@test size(LARLIB.sphere(2.5,pi/3,pi/5)()[1],2)==703
	@test length(LARLIB.sphere(2.5,pi/3,pi/5)()[2])==1296
end



@testset "toroidal" begin
	@test BoxCalculation(LARLIB.toroidal(1,3,2*pi,2*pi)()[1])==128
	@test BoxCalculation(LARLIB.toroidal(2,3,2*pi,2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.toroidal(1.3,4.6,pi/4,pi/7)()[1],2)==925
	@test length(LARLIB.toroidal(1.3,4.6,pi/4,pi/7)()[2])==1728
end



@testset "crown" begin
	@test BoxCalculation(LARLIB.crown(1., 3., 2*pi)()[1])==128
	@test BoxCalculation(LARLIB.crown(2., 3., 2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.crown(1.5, 5.6, pi/8)()[1],2)==481
	@test length(LARLIB.crown(1.5, 5.6, pi/8)()[2])==864
end




@testset "ball" begin
	@test BoxCalculation(LARLIB.ball(1, pi, 2*pi)()[1])==8
	@test BoxCalculation(LARLIB.ball(6, pi, pi)()[1])==864
	@test BoxCalculation(LARLIB.ball(3, pi, 2*pi)()[1])==6^3
	@test size(LARLIB.ball(2.6, pi/5, pi/9)()[1],2)==2813
	@test length(LARLIB.ball(2.6, pi/5, pi/9)()[2])==2592
end



@testset "rod" begin
	@test BoxCalculation(LARLIB.rod(1,5,2*pi)()[1])==20
	@test BoxCalculation(LARLIB.rod(2,2,pi)()[1])==16
	@test BoxCalculation(LARLIB.rod(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(LARLIB.rod(3.7,8.9,pi/9)()[1],2)==74
	@test length(LARLIB.rod(3.7,8.9,pi/9)()[2])==1
end



@testset "hollowCyl" begin
	@test BoxCalculation(LARLIB.hollowCyl(0,1.,5,2*pi)()[1])==20
	@test BoxCalculation(LARLIB.hollowCyl(1,2.,4,pi)()[1])==32
	@test BoxCalculation(LARLIB.hollowCyl(1,4,3,pi/2)()[1])==48
	#((radius*2)^2)*height
	@test size(LARLIB.hollowCyl(3,4.,7.8,pi/5)()[1],2)==148
	@test length(LARLIB.hollowCyl(3,4.,7.8,pi/5)()[2])==36
end



@testset "hollowBall" begin
	@test BoxCalculation(LARLIB.hollowBall(1,2,pi,2*pi)([36,36,1])[1]) ==64
	@test BoxCalculation(LARLIB.hollowBall(1,6,pi,pi)([36,36,1])[1])== 864
	@test BoxCalculation(LARLIB.hollowBall(2,4,pi,2*pi)([36,36,1])[1])==8^3
	#(radius*2)^3
	@test size(LARLIB.hollowBall(1.5,6.7,pi/3,2*pi/3)()[1],2)==3700
	@test length(LARLIB.hollowBall(1.5,6.7,pi/3,2*pi/3)()[2])==2592
end



@testset "torus" begin
	@test BoxCalculation(LARLIB.torus(1.,2.,.5,2*pi,2*pi)()[1])==147
	@test BoxCalculation(LARLIB.torus(2,3,.5,2*pi,2*pi)()[1])==605.0
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.torus(5.2,7,.5,pi/3,pi/4)()[1],2)==4625
	@test length(LARLIB.torus(5.2,7,.5,pi/3,pi/4)()[2])==3456
end







