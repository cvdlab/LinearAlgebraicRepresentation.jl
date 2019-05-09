using LinearAlgebraicRepresentation
using Plasm, SparseArrays, LinearAlgebra
Lar = LinearAlgebraicRepresentation
using Rebugger, Triangle
using JuliaInterpreter

# source data from minimal test example
v = [0.0 1.0 0.0 0.0 0.5 0.25 0.25 0.25 1.25 0.25 0.25; 0.0 0.0 1.0 0.0 0.25 0.5 0.25 0.25 0.25 1.25 0.25; 0.0 0.0 0.0 1.0 0.25 0.25 0.5 0.25 0.25 0.25 1.25]
ev = Array{Int64,1}[[1, 2], [1, 3], [2, 3], [1, 4], [2, 4], [3, 4], [5, 6], [5, 7], [6, 7], [5, 8], [5, 9], [6, 8], [6, 10], [9, 10], [7, 8], [7, 11], [9, 11], [10, 11]]
cscEV = ([1, 2, 4, 1, 3, 5, 2, 3, 6, 4, 5, 6, 7, 8, 10, 11, 7, 9, 12, 13, 8, 9, 15, 16, 10, 12, 15, 11, 14, 17, 13, 14, 18, 16, 17, 18], [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11], Int8[-1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1])
fe = Array{Int64,1}[[1, 2, 3], [1, 4, 5], [2, 4, 6], [3, 5, 6, 7, 8, 9], [7, 8, 9], [7, 10, 12], [7, 11, 13, 14], [8, 10, 15], [8, 11, 16, 17], [9, 12, 15], [9, 13, 16, 18], [14, 17, 18]]
cscFE = ([1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4, 4, 5, 6, 7, 4, 5, 8, 9, 4, 5, 10, 11, 6, 8, 7, 9, 6, 10, 7, 11, 7, 12, 8, 10, 9, 11, 9, 12, 11, 12], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18], Int8[1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1, 1])
fv = Array{Int64,1}[[2, 3, 1], [4, 2, 1], [4, 3, 1], [7, 4, 2, 3, 5, 6], [7, 5, 6], [5, 8, 6], [9, 10, 5, 6], [7, 5, 8], [7, 9, 11, 5], [7, 8, 6], [7, 10, 11, 6], [9, 10, 11]]

Plasm.view(Plasm.numbering(0.25)((v,[[[k] for k=1:size(v,2)],ev,fv])))

#V, copEV, copFE, copCF = Lar.Arrangement.spatial_arrangement( W::Lar.Points, cop_EW::Lar.ChainOp, cop_FE::Lar.ChainOp)

copEV = SparseArrays.sparse(cscEV...);
copFE = SparseArrays.sparse(cscFE...);
V = convert(Lar.Points, v');

#breakpoint(Triangle.constrained_triangulation)

function minimal_3cycles(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)

    triangulated_faces = Array{Any, 1}(undef, FE.m)

    function face_angle(e::Int, f::Int)
        if !isassigned(triangulated_faces, f)
            vs_idxs = Array{Int64, 1}()
            edges_idxs = FE[f, :].nzind
            edge_num = length(edges_idxs)
            edges = zeros(Int64, edge_num, 2)

            for (i, ee) in enumerate(edges_idxs)
                edge = EV[ee, :].nzind
                edges[i, :] = edge
                vs_idxs = union(vs_idxs, edge)
            end

            vs = V[vs_idxs, :]

            v1 = normalize(vs[2, :] - vs[1, :])
            v2 = [0 0 0]		# added for debug
            v3 = [0 0 0]
            err = 1e-8
            i = 3
            while -err < norm(v3) < err
                v2 = normalize(vs[i, :] - vs[1, :])
                v3 = cross(v1, v2)
                i = i + 1
            end

            M = reshape([v1; v2; v3], 3, 3)

            vs = vs*M

            triangulated_faces[f] = Triangle.constrained_triangulation(
                Array{Float64,2}(vs), vs_idxs, edges, fill(true, edge_num))

        end
        edge_vs = EV[e, :].nzind

        t = findfirst(x->edge_vs[1] in x && edge_vs[2] in x, triangulated_faces[f])



        v1 = normalize(V[edge_vs[2], :] - V[edge_vs[1], :])

        if abs(v1[1]) > abs(v1[2])
            invlen = 1. / sqrt(v1[1]*v1[1] + v1[3]*v1[3])
            v2 = [-v1[3]*invlen, 0, v1[1]*invlen]
        else
            invlen = 1. / sqrt(v1[2]*v1[2] + v1[3]*v1[3])
            v2 = [0, -v1[3]*invlen, v1[2]*invlen]
        end

        v3 = cross(v1, v2)

        M = reshape([v1; v2; v3], 3, 3)

        triangle = triangulated_faces[f][t]
        third_v = setdiff(triangle, edge_vs)[1]
        vs = V[[edge_vs..., third_v], :]*M

        v = vs[3, :] - vs[1, :]
        angle = atan(v[2], v[3])
        return angle
    end

    #EF = FE'
    EF = convert(Lar.ChainOp, LinearAlgebra.transpose(FE))
@show face_angle;
@show V;
@show findnz(EF);
	FC = Lar.Arrangement.minimal_cycles(face_angle, true)(V, EF)

	#FC'
    return -convert(Lar.ChainOp, LinearAlgebra.transpose(FC))
end

# Juno.@enter debug_function()
function debug_function()
	minimal_3cycles(V::Lar.Points, copEV::Lar.ChainOp, copFE::Lar.ChainOp)

end
# minimal_3cycles(V::Lar.Points, copEV::Lar.ChainOp, copFE::Lar.ChainOp)
