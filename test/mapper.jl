
using Base.Test
using LARLIB.Test

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
	@test BoxCalculation(LARLIB.circle(radius=2.;angle=2*pi)()[1])==16
	@test BoxCalculation(LARLIB.circle(radius=3;angle=2*pi)()[1])==36
	@test BoxCalculation(LARLIB.circle(radius=5,angle=pi/2)()[1])==25
	@test size(LARLIB.circle(radius=3,angle=2*pi)(60)[1],2)==60
	@test length(LARLIB.circle(radius=3,angle=2*pi)(60)[2])==60
end



@testset "helix" begin
	@test BoxCalculation(LARLIB.helix()()[1])==8
	@test BoxCalculation(LARLIB.helix(radius=2, pitch=2, nturns=2)()[1])==64
	@test BoxCalculation(LARLIB.helix(radius=1, pitch=2, nturns=5)()[1])==40
	@test BoxCalculation(LARLIB.helix(radius=1, pitch=1, nturns=3)()[1])==12
	#=
	radius=1, pitch=1, nturns=2, is contained
	in a box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LARLIB.helix(radius=5, pitch=7, nturns=9)()[1],2)==325
	@test length(LARLIB.helix(radius=5, pitch=7, nturns=9)()[2])==324
end



@testset "disk" begin
	@test BoxCalculation(LARLIB.disk(radius=1., angle=2*pi)([36, 1])[1])==4
	@test BoxCalculation(LARLIB.disk(radius=2., angle=2*pi)()[1])==16
	@test BoxCalculation(LARLIB.disk(radius=2., angle=pi)()[1])==8
	@test BoxCalculation(LARLIB.disk(radius=2., angle=pi/2)()[1])==4
	@test size(LARLIB.disk(radius=10, angle=pi/7)()[1],2)==75
	@test length(LARLIB.disk(radius=10, angle=pi/7)()[2])==108
end



@testset "helicoid" begin
	@test BoxCalculation(LARLIB.helicoid()()[1])==8
	@test BoxCalculation(LARLIB.helicoid(R=2,r=1,pitch=2,nturns=2)()[1])==64
	@test BoxCalculation(LARLIB.helicoid(R=1,r=.5,pitch=2,nturns=5)()[1])==40
	@test BoxCalculation(LARLIB.helicoid(R=1,r=.3,pitch=1,nturns=3)()[1])==12
	#=
	radius=1, pitch=1, nturns=2, is contained in a
	box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(LARLIB.helicoid(R=1.3,r=.7,pitch=1,nturns=3)()[1],2)==327
	@test length(LARLIB.helicoid(R=1.3,r=.7,pitch=1,nturns=3)()[2])==432
end


@testset "ring" begin
	@test BoxCalculation(LARLIB.ring(r=1.,R=3.,angle=2*pi)()[1])==36
	@test BoxCalculation(LARLIB.ring(r=1.,R=2,angle=pi)()[1])==8
	@test BoxCalculation(LARLIB.ring(r=1.,R=5,angle=pi/2)()[1])==25
	#(radius*2)^2
	@test size(LARLIB.ring(r=1,R=5,angle=pi/2)()[1],2)==74
	@test length(LARLIB.ring(r=5,R=10,angle=pi/6)()[2])==36
end



@testset "cylinder" begin
	@test BoxCalculation(LARLIB.cylinder(radius=1,height=5,angle=2*pi)()[1])==20
	@test BoxCalculation(LARLIB.cylinder(radius=2,height=2,angle=pi)()[1])==16
	@test BoxCalculation(LARLIB.cylinder(radius=1,height=4,angle=pi)()[1])==8
	#((radius*2)^2)*height
	@test size(LARLIB.cylinder(radius=3.4,height=20,angle=pi/7)()[1],2)==74
	@test length(LARLIB.cylinder(radius=3.4,height=20,angle=pi/7)()[2])==36
end



@testset "sphere" begin
	@test BoxCalculation(LARLIB.sphere(radius=2,angle1=pi,angle2=2*pi)()[1])==64
	@test BoxCalculation(LARLIB.sphere(radius=6,angle1=pi,angle2=pi)()[1])==864
	@test BoxCalculation(LARLIB.sphere(radius=4,angle1=pi,angle2=2*pi)()[1])==8^3
	@test size(LARLIB.sphere(radius=2.5,angle1=pi/3,angle2=pi/5)()[1],2)==703
	@test length(LARLIB.sphere(radius=2.5,angle1=pi/3,angle2=pi/5)()[2])==1296
end



@testset "toroidal" begin
	@test BoxCalculation(LARLIB.toroidal(r=1,R=3,angle1=2*pi,angle2=2*pi)()[1])==128
	@test BoxCalculation(LARLIB.toroidal(r=2,R=3,angle1=2*pi,angle2=2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.toroidal(r=1.3,R=4.6,angle1=pi/4,angle2=pi/7)()[1],2)==925
	@test length(LARLIB.toroidal(r=1.3,R=4.6,angle1=pi/4,angle2=pi/7)()[2])==1728
end



@testset "crown" begin
	@test BoxCalculation(LARLIB.crown(r=1., R=3., angle=2*pi)()[1])==128
	@test BoxCalculation(LARLIB.crown(r=2., R=3., angle=2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.crown(r=1.5, R=5.6, angle=pi/8)()[1],2)==481
	@test length(LARLIB.crown(r=1.5, R=5.6, angle=pi/8)()[2])==864
end




@testset "ball" begin
	@test BoxCalculation(LARLIB.ball(radius=1, angle1=pi, angle2=2*pi)()[1])==8
	@test BoxCalculation(LARLIB.ball(radius=6, angle1=pi, angle2=pi)()[1])==864
	@test BoxCalculation(LARLIB.ball(radius=3, angle1=pi, angle2=2*pi)()[1])==6^3
	@test size(LARLIB.ball(radius=2.6, angle1=pi/5, angle2=pi/9)()[1],2)==2813
	@test length(LARLIB.ball(radius=2.6, angle1=pi/5, angle2=pi/9)()[2])==2592
end



@testset "rod" begin
	@test BoxCalculation(LARLIB.rod(radius=1,height=5,angle=2*pi)()[1])==12
	@test BoxCalculation(LARLIB.rod(radius=2,height=2,angle=pi)()[1])==12
	@test BoxCalculation(LARLIB.rod(radius=1,height=4,angle=pi)()[1])==12
	#((radius*2)^2)*height
	@test size(LARLIB.rod(radius=3.7,height=8.9,angle=pi/9)()[1],2)==72
	@test length(LARLIB.rod(radius=3.7,height=8.9,angle=pi/9)()[2])==1
end



@testset "hollowCyl" begin
	@test BoxCalculation(LARLIB.hollowCyl(r=0,R=1.,height=5,angle=2*pi)()[1])==20
	@test BoxCalculation(LARLIB.hollowCyl(r=1,R=2.,height=4,angle=pi)()[1])==32
	@test BoxCalculation(LARLIB.hollowCyl(r=1,R=4,height=3,angle=pi/2)()[1])==48
	#((radius*2)^2)*height
	@test size(LARLIB.hollowCyl(r=3,R=4.,height=7.8,angle=pi/5)()[1],2)==148
	@test length(LARLIB.hollowCyl(r=3,R=4.,height=7.8,angle=pi/5)()[2])==36
end



@testset "hollowBall" begin
	@test BoxCalculation(LARLIB.hollowBall(r=1,R=2,angle1=pi,angle2=2*pi)([36,36,1])[1]) ==64
	@test BoxCalculation(LARLIB.hollowBall(r=1,R=6,angle1=pi,angle2=pi)([36,36,1])[1])== 864
	@test BoxCalculation(LARLIB.hollowBall(r=2,R=4,angle1=pi,angle2=2*pi)([36,36,1])[1])==8^3
	#(radius*2)^3
	@test size(LARLIB.hollowBall(r=1.5,R=6.7,angle1=pi/3,angle2=2*pi/3)()[1],2)==3700
	@test length(LARLIB.hollowBall(r=1.5,R=6.7,angle1=pi/3,angle2=2*pi/3)()[2])==2592
end



@testset "torus" begin
	@test BoxCalculation(LARLIB.torus(r=1.,R=2.,h=.5,angle1=2*pi,angle2=2*pi)()[1])==147
	@test BoxCalculation(LARLIB.torus(r=2,R=3,h=.5,angle1=2*pi,angle2=2*pi)()[1])==605.0
	#(((R+r)*2)^2)*(r*2)
	@test size(LARLIB.torus(r=5.2,R=7,h=.5,angle1=pi/3,angle2=pi/4)()[1],2)==4625
	@test length(LARLIB.torus(r=5.2,R=7,h=.5,angle1=pi/3,angle2=pi/4)()[2])==3456
end







