
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
	@test BoxCalculation(circle()()[1])==4
	@test BoxCalculation(circle(2,2*pi)()[1])==16
	@test BoxCalculation(circle(3,pi)()[1])==18
	@test BoxCalculation(circle(5,pi/2)()[1])==25
	#(radius*2)^2
	@test size(circle(3,2*pi)(60)[1],2)==60
	@test length(circle(3,2*pi)(60)[2])==60
	#from python
end



@testset "helix" begin
	@test BoxCalculation(helix()()[1])==8
	@test BoxCalculation(helix(2,2,2)()[1])==64
	@test BoxCalculation(helix(1,2,5)()[1])==40
	@test BoxCalculation(helix(1,1,3)()[1])==12
	#=
	radius=1, pitch=1, nturns=2, is contained
	in a box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(helix(5,7,9)()[1],2)==325
	@test length(helix(5,7,9)()[2])==324
end



@testset "disk" begin
	@test BoxCalculation(disk()()[1])==4
	@test BoxCalculation(disk(2,2*pi)()[1])==16
	@test BoxCalculation(disk(3,pi)()[1])==18
	@test BoxCalculation(disk(5,pi/2)()[1])==25
	#(radius*2)^2
	@test size(disk(10,pi/7)()[1],2)==75
	@test length(disk(10,pi/7)()[2])==108
end



@testset "helicoid" begin
	@test BoxCalculation(helicoid()()[1])==8
	@test BoxCalculation(helicoid(2,1,2,2)()[1])==64
	@test BoxCalculation(helicoid(1,0.5,2,5)()[1])==40
	@test BoxCalculation(helicoid(1,0.3,1,3)()[1])==12
	#=
	radius=1, pitch=1, nturns=2, is contained in a
	box that has volume (radius*2)^2*pitch*nturns
	=#
	@test size(helicoid(1.3,0.7,1,3)()[1],2)==327
	@test length(helicoid(1.3,0.7,1,3)()[2])==432
end


@testset "ring" begin
	@test BoxCalculation(ring(1,3,2*pi)()[1])==36
	@test BoxCalculation(ring(1,2,pi)()[1])==8
	@test BoxCalculation(ring(2,5,pi/2)()[1])==25
	#(radius*2)^2
	@test size(ring(5,10,pi/6)()[1],2)==74
	@test length(ring(5,10,pi/6)()[2])==36
end



@testset "cylinder" begin
	@test BoxCalculation(cylinder(1,5,2*pi)()[1])==20
	@test BoxCalculation(cylinder(2,2,pi)()[1])==16
	@test BoxCalculation(cylinder(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(cylinder(3.4,20,pi/7)()[1],2)==74
	@test length(cylinder(3.4,20,pi/7)()[2])==36
end



@testset "sphere" begin
	@test BoxCalculation(sphere(2,pi,2*pi)()[1])==64
	@test BoxCalculation(sphere(6,pi,pi)()[1])==864
	@test BoxCalculation(sphere(4,pi,2*pi)()[1])==8^3
	#(radius*2)^3
	@test size(sphere(2.5,pi/3,pi/5)()[1],2)==703
	@test length(sphere(2.5,pi/3,pi/5)()[2])==1296
end



@testset "toroidal" begin
	@test BoxCalculation(toroidal(1,3,2*pi,2*pi)()[1])==128
	@test BoxCalculation(toroidal(2,3,2*pi,2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(toroidal(1.3,4.6,pi/4,pi/7)()[1],2)==925
	@test length(toroidal(1.3,4.6,pi/4,pi/7)()[2])==1728
end



@testset "crown" begin
	@test BoxCalculation(crown(1,3,2*pi)()[1])==128
	@test BoxCalculation(crown(2,3,2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(crown(1.5,5.6,pi/8)()[1],2)==481
	@test length(crown(1.5,5.6,pi/8)()[2])==864
end




@testset "ball" begin
	@test BoxCalculation(ball(1,pi,2*pi)()[1])==8
	@test BoxCalculation(ball(6,pi,pi)()[1])==864
	@test BoxCalculation(ball(3,pi,2*pi)()[1])==6^3
	@test size(ball(2.6,pi/5,pi/9)()[1],2)==2813
	@test length(ball(2.6,pi/5,pi/9)()[2])==2592
end



@testset "rod" begin
	@test BoxCalculation(rod(1,5,2*pi)()[1])==20
	@test BoxCalculation(rod(2,2,pi)()[1])==16
	@test BoxCalculation(rod(1,4,pi)()[1])==8
	#((radius*2)^2)*height
	@test size(rod(3.7,8.9,pi/9)()[1],2)==74
	@test length(rod(3.7,8.9,pi/9)()[2])==1
end



@testset "hollowCyl" begin
	@test BoxCalculation(hollowCyl(0,1,5,2*pi)()[1])==20
	@test BoxCalculation(hollowCyl(1,2,4,pi)()[1])==32
	@test BoxCalculation(hollowCyl(1,4,3,pi/2)()[1])==48
	#((radius*2)^2)*height
	@test size(hollowCyl(3.4,7.8,pi/5)()[1],2)==144
	@test length(hollowCyl(3.4,7.8,pi/5)()[2])==36
end



@testset "hollowBall" begin
	@test BoxCalculation(hollowBall(0,2,pi,2*pi)([36,36,1])[1])==64
	@test BoxCalculation(hollowBall(1,6,pi,pi)([36,36,1])[1])==864
	@test BoxCalculation(hollowBall(2,4,pi,2*pi)([36,36,1])[1])==8^3
	#(radius*2)^3
	@test size(hollowBall(1.5,6.7,pi/3,2*pi/3)()[1],2)==3700
	@test length(hollowBall(1.5,6.7,pi/3,2*pi/3)()[2])==2592
end



@testset "torus" begin
	@test BoxCalculation(torus(1,3,2*pi,2*pi)()[1])==128
	@test BoxCalculation(torus(2,3,2*pi,2*pi)()[1])==400
	#(((R+r)*2)^2)*(r*2)
	@test size(torus(5.2,7,pi/3,pi/4)()[1],2)==3737
	@test length(torus(5.2,7,pi/3,pi/4)()[2])==3456
end



@testset "pizza" begin
	@test BoxCalculation(pizza(1,2,pi/2)()[1])==18
	@test BoxCalculation(pizza(0.5,2,pi)()[1])==12.5
	#(((R+r)*2)^2)*r if angle=pi
	@test size(pizza(0.03,5.6,pi/9)()[1],2)==927
	@test length(pizza(0.03,5.6,pi/9)()[2])==1
end






