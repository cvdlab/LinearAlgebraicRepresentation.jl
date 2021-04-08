#
#1. Create a regular $N^3$ grid T of random tetrahedra:
#	(a) Generate four vertices for each tetrahedron randomly in a half-open fixed-size cube.
#	(b) Let these cubes form a regular $N^3$ grid.
#2. Create a regular $(N − 1)^3$ grid C of such cubes.
#3. Align T and C such that the grid nodes of C are at the centers of the grid cells of T .
#4. Measure time for T ∪ C.
#

# TODO: DEBUG !!!
#----------------

using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation
using ViewerGL
GL = ViewerGL


function t1t2(i,j,k)
	v1 = rand(3,4).+ [2i,2j,2k]
	v2 = rand(3,4).+ [2i,2j,2k]
	ev = [[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
	fv = [[2,3,4],[1,3,4],[1,2,4],[1,2,3]]
	t1 = [v1,fv,ev] 
	t2 = [v2,fv,ev]
	return t1,t2
end

N = 2  # N ≥ 2
tets = []
for i=0:N-1, j=0:N-1, k=0.5:0.5
	push!(tets, t1t2(i,j,k)...)
end
V,FV,EV = Lar.struct2lar(Lar.Struct(tets))
GL.VIEW([ GL.GLGrid(V,FV, GL.COLORS[1],1), GL.GLFrame2 ]);

W,(_,EW,FW,_) = Lar.cuboidGrid([1,1,1],true);
U = (W .* [2(N),2(N),2]) .+ [-.5,-.5,-.5]
GL.VIEW([ GL.GLGrid(V,FV,GL.COLORS[1],0.5), GL.GLGrid(U,FW,GL.COLORS[2],0.5), GL.GLFrame2 ]);

V,FV,EV = Lar.struct2lar(Lar.Struct([Lar.Struct(tets), Lar.Struct([(U,FW,EW)]) ]))
GL.VIEW([ GL.GLGrid(V,EV, GL.COLORS[1],0.5), GL.GLFrame2 ]);

	cop_EV = convert(Lar.ChainOp, Lar.coboundary_0(EV::Lar.Cells));
	cop_FE = Lar.coboundary_1(V, FV::Lar.Cells, EV::Lar.Cells); 
	W = convert(Lar.Points, V');
	V, copEV, copFE, copCF = Lar.space_arrangement( W, cop_EV, cop_FE);
	W = convert(Lar.Points, V');
	V,CVs,FVs,EVs = Lar.pols2tria(W, copEV, copFE, copCF);
	
GL.VIEW( GL.GLExplode(V,CVs[3:end],1.5,1.5,1.5,99,0.5) );
GL.VIEW( GL.GLExplode(V,CVs[3:end],1,1,1,99,0.5) );
GL.VIEW( GL.GLExplode(V,CVs[2:2],5,5,5,99,0.5) );
GL.VIEW( GL.GLExplode(V,CVs[12:12],5,5,5,99,0.5) );
