using Triangulate

"""
    constrained_triangulation2D(V::Lar.Points, EV::Lar.Cells) -> Lar.Cells
"""
function constrained_triangulation2D(V::Lar.Points, EV::Lar.Cells)
	triin = Triangulate.TriangulateIO()
	triin.pointlist = V
	triin.segmentlist = hcat(EV...)
	(triout, vorout) = triangulate("pQ", triin)
	trias = Array{Int64,1}[c[:] for c in eachcol(triout.trianglelist)]
	return trias
end


"""
 	LinearAlgebraicRepresentation.triangulate2d(V::Lar.Points, EV::Lar.Cells)

@overwrite Lar.triangulate2d.
"""
function LinearAlgebraicRepresentation.triangulate2d(V::Lar.Points, EV::Lar.Cells)
   	 # data for Constrained Delaunay Triangulation (CDT)
   	 points = convert(Array{Float64,2}, V')
	 # points_map = Array{Int64,1}(collect(1:1:size(points)[1]))
   	 # edges_list = convert(Array{Int64,2}, hcat(EV...)')
   	 # edge_boundary = [true for k=1:size(edges_list,1)] ## dead code !!

	trias = constrained_triangulation2D(V::Lar.Points, EV::Lar.Cells)

 	#Triangle.constrained_triangulation(points,points_map,edges_list)
	innertriangles = Array{Int64,1}[]
	for (u,v,w) in trias
		point = (points[u,:]+points[v,:]+points[w,:])./3
		copEV = Lar.lar2cop(EV)
		inner = Lar.point_in_face(point, points::Lar.Points, copEV::Lar.ChainOp)
		if inner
			push!(innertriangles,[u,v,w])
		end
	end
    return innertriangles
end
