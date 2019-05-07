function minimal_cycles(angles_fn::Function, verbose=false)

    function _minimal_cycles(V::Lar.Points,
    ld_bounds::Lar.ChainOp)

@show V
@show ld_bounds

        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int8, ld_cellsnum)
        dir_marks = zeros(Int8, ld_cellsnum)
        d_bounds = spzeros(Int8, ld_cellsnum, 0)

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
@show angles
        function nextprev(lld::Int64, ld::Int64, norp)
            as = angles[lld]
            #ne = findfirst(as, ld)  (findfirst(isequal(v), A), 0)[1]
            ne = (findfirst(isequal(ld), as), 0)[1]
            while true
                ne += norp
                if ne > length(as)
                    ne = 1
                elseif ne < 1
                    ne = length(as)
                end

                if count_marks[as[ne]] < 2
                    break
                end
            end
            as[ne]
        end


        while (sigma = get_seed_cell()) > 0
@show sigma

            if verbose
                print(Int(floor(50 * sum(count_marks) / ld_cellsnum)), "%\r")
            end

            c_ld = spzeros(Int8, ld_cellsnum)
            if count_marks[sigma] == 0
                c_ld[sigma] = 1
            else
                c_ld[sigma] = -dir_marks[sigma]
            end
            c_lld = ld_bounds*c_ld
            while c_lld.nzind != []
				#corolla = spzeros(Int8, ld_cellsnum)
				corolla = zeros(Int64, ld_cellsnum)
                for tau in c_lld.nzind
                    b_ld = ld_bounds[tau, :]
                    pivot = intersect(c_ld.nzind, b_ld.nzind)[1]
                    adj = nextprev(tau, pivot, sign(-c_lld[tau]))
                    corolla[adj] = c_ld[pivot]
                    if b_ld[adj] == b_ld[pivot]
                        corolla[adj] *= -1
                    end
@show corolla
                end
                c_ld += corolla
                c_lld = ld_bounds*c_ld
@show c_lld
            end
            map(s->count_marks[s] += 1, c_ld.nzind)
            map(s->dir_marks[s] = c_ld[s], c_ld.nzind)
            d_bounds = [d_bounds c_ld]
@show d_bounds
        end
        return d_bounds

    end

    return _minimal_cycles
end
