include("./utilities.jl")

function minimal_2cycles(V::Verts, ev::Cells)

    function edge_angle(V::Verts, v::Int, edge::Cell)
        v2 = setdiff(edge.nzind, [v])[1]
        x, y = V[v2, :] - V[v, :]
        return atan2(y, x)
    end

    for i in 1:ev.m
        j = ev[i,:].nzind[1]
        ev[i, j] = -1
    end
    VE = ev'

    EF = minimal_cycles(edge_angle)(V, VE)

    return EF'
end


function minimal_cycles(angles_fn::Function)

    function _minimal_cycles(V::Verts, ld_bounds::Cells)
        lld_cellsnum, ld_cellsnum = size(ld_bounds)
        count_marks = zeros(Int8, ld_cellsnum)
        dir_marks = zeros(Int8, ld_cellsnum)
        d_bounds = spzeros(Int8, ld_cellsnum, 0)
        
        angles = Array{Array{Int64, 1}, 1}(lld_cellsnum)
        
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
            as = [(ld, angles_fn(V, lld, ld_bounds[:, ld])) 
                for ld in ld_bounds[lld, :].nzind]
            sort!(as, lt=(a,b)->a[2]<b[2])
            as = map(a->a[1], as)
            angles[lld] = as
        end
        function nextprev(lld::Int64, ld::Int64, norp)
            as = angles[lld]
            ne = findfirst(as, ld)
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
            c_ld = spzeros(Int8, ld_cellsnum)
            if count_marks[sigma] == 0
                c_ld[sigma] = 1
            else
                c_ld[sigma] = -dir_marks[sigma]
            end
            c_lld = ld_bounds*c_ld
            while c_lld.nzind != []
                corolla = spzeros(Int8, ld_cellsnum)
                for tau in c_lld.nzind
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
            map(s->dir_marks[s] = c_ld[s], c_ld.nzind)
            d_bounds = [d_bounds c_ld]
            
        end
        
        return d_bounds
        
    end

    return _minimal_cycles
end
