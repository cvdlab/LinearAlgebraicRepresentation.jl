#=
1. Create set of N triangles T such that
	(a) they all meet in a common vertex v.
	(b) the edges E opposite to v are parallel to each other with uniform distance between 
	pairs of triangles. 
	(c) the length of each edge e in E equals the largest distance between to pairs of edges 		in E.
2. Create a copy T′ of T.
3. Rotate T′ around the axis through v and the center of the edges E by 90 degrees.
4. Measure time for T ∪ T′.
=#

using ViewerGL; GL = ViewerGL
using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation;
