using LinearAlgebraicRepresentation
Lar = LinearAlgebraicRepresentation

function minimal_2cycles(V::Lar.Points, EV::Lar.ChainOp)

    function edge_angle(v::Int, e::Int)
        edge = EV[e, :]
        v2 = setdiff(edge.nzind, [v])[1]
        x, y = V[v2, :] - V[v, :]
        return atan(y, x)
    end

    for i in 1:EV.m
        j = min(EV[i,:].nzind...)
        EV[i, j] = -1
    end
    VE = convert(Lar.ChainOp, SparseArrays.transpose(EV))
    EF = Lar.Arrangement.minimal_cycles(edge_angle)(V, VE)

    return convert(Lar.ChainOp, SparseArrays.transpose(EF))
end



function minimal_3cycles(V::Lar.Points, EV::Lar.ChainOp, FE::Lar.ChainOp)
		@show V;
		@show findnz(EV);
		@show findnz(FE);

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

            #vs = V[vs_idxs, :]
			fv,edges = Lar.vcycle(EV, FE, f)

			vs = V[fv, :]


            v1 = LinearAlgebra.normalize(vs[2, :] - vs[1, :])
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

            #vs = vs*M
			vs = (vs*M)[:, 1:2]

            # triangulated_faces[f] = Triangle.constrained_triangulation(
            #     Array{Float64,2}(vs), vs_idxs, edges, fill(true, edge_num))
			v = convert(Lar.Points, vs'[1:2,:])
			vmap = Dict(zip(fv,1:length(fv))) # vertex map
			mapv = Dict(zip(1:length(fv),fv)) # inverse vertex map
			trias = Lar.triangulate2d(v,edges)
			triangulated_faces[f] = [[mapv[v] for v in tria] for tria in trias]
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
		println(">>>>>>>> going to loop?")
		FC = Lar.Arrangement.minimal_cycles(face_angle, true)(V, EF)  # , EV)

		#FC'
    return -convert(Lar.ChainOp, LinearAlgebra.transpose(FC))
end


function minimal_cycles(angles_fn::Function, verbose=true)

    function _minimal_cycles(V::Lar.Points,
    ld_bounds::Lar.ChainOp)  # , EV)

        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int64, ld_cellsnum)
        dir_marks = zeros(Int64, ld_cellsnum)
        d_bounds = spzeros(Int64, ld_cellsnum, 0)

        angles = Array{Array{Int64, 1}, 1}(undef, lld_cellsnum)

        function get_seed_cell()
            s = -1
            for i in 1:ld_cellsnum
                if count_marks[i] == 0
                    return i
                elseif count_marks[i] == 1 && s < 0
                    s = i
                end
            end
            return s
        end
        for lld in 1:lld_cellsnum
            as = []
            for ld in ld_bounds[lld, :].nzind
                push!(as, (ld, angles_fn(lld, ld)))
            end
            sort!(as, lt=(a,b)->a[2]<b[2])
            as = map(a->a[1], as)
            angles[lld] = as
        end
        function nextprev(lld::Int64, ld::Int64, norp)
            as = angles[lld]
            #ne = findfirst(as, ld)  (findfirst(isequal(v), A), 0)[1]
            ne = (findfirst(isequal(ld), as), 0)[1]
            while true
                ne += norp # next or previous
                if ne > length(as)
                    ne = 1
                elseif ne < 1
                    ne = length(as)
                end

                if count_marks[as[ne]] < 2
                    break
                end
            end
            for k=1:length(count_marks)
            	if count_marks[k]>2  error("TGW is looping") end
            end

            as[ne]
        end

        while (sigma = get_seed_cell()) > 0
						if verbose
                println(Int(floor(50 * sum(count_marks) / ld_cellsnum)), "%\r") # <<<<<<<<<<<<<<<<<<<
            end

            c_ld = spzeros(Int8, ld_cellsnum)
            if count_marks[sigma] == 0
                c_ld[sigma] = 1
            else
                c_ld[sigma] = -dir_marks[sigma]
            end
            c_lld = ld_bounds*c_ld
            while c_lld.nzind != []
            		println(">>>>>> eccomi!")
								corolla = spzeros(Int64, ld_cellsnum)
								#corolla = zeros(Int64, ld_cellsnum)
								
                for tau in c_lld.nzind # when looping, loops here !!
                @show tau,c_lld.nzind
                    b_ld = ld_bounds[tau, :]
                    pivot = intersect(c_ld.nzind, b_ld.nzind)[1]
                    adj = nextprev(tau, pivot, sign(-c_lld[tau]))
                    corolla[adj] = c_ld[pivot]
                    if b_ld[adj] == b_ld[pivot]
                        corolla[adj] *= -1
                    end
                end
                
                c_ld += corolla
                c_lld = ld_bounds*c_ld
            end
            map(s->count_marks[s] += 1, c_ld.nzind)
            @show count_marks
            for k=1:length(count_marks)
            	if count_marks[k]>2  error("TGW is looping") end
            end
            map(s->dir_marks[s] = c_ld[s], c_ld.nzind)
            @show c_ld; println("")
            d_bounds = [d_bounds c_ld]
            
#           V = convert(Lar.Points, V');
#						V,CVs,FVs,EVs = Lar.pols2tria(V, EV, ld_bounds, d_bounds) # partial assembly
#						GL.VIEW(GL.GLExplode(V,FVs,1.1,1.1,1.1,99,1));
#						GL.VIEW(GL.GLExplode(V,EVs,1.5,1.5,1.5,99,1));
#						GL.VIEW(GL.GLExplode(V,CVs,1,1,1,99,1));

        end
        return d_bounds

    end

    return _minimal_cycles
end
