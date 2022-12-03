using LinearAlgebraicRepresentation
using SparseArrays
Lar = LinearAlgebraicRepresentation
using Debugger


V = [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0; 0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5; 0.25 0.25 0.25; 1.25 0.25 0.25; 0.25 1.25 0.25; 0.25 0.25 1.25]
ld_bounds = ([1, 2, 3, 1, 4, 5, 2, 4, 6, 3, 5, 6, 7, 8, 9, 7, 8, 9, 7, 10, 12, 7, 11, 13, 14, 8, 10, 15, 8, 11, 16, 17, 9, 12, 15, 9, 13, 16, 18, 14, 17, 18], [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12], Int8[1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1])
angles = Array{Int64,1}[[2, 1], [1, 3], [4, 1], [3, 2], [2, 4], [4, 3], [7, 4, 5, 6], [9, 8, 4, 5], [11, 4, 5, 10], [8, 6], [9, 7], [10, 6], [7, 11], [7, 12], [8, 10], [11, 9], [9, 12], [12, 11]]

sigma =1

function get_seed_cell(ld_cellsnum, count_marks)
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

function _minimal_cycles(V::Lar.Points,ld_bounds::Lar.ChainOp)
#@show V
#@show ld_bounds
    lld_cellsnum, ld_cellsnum = size(ld_bounds)
    count_marks = zeros(Int8, ld_cellsnum)
    dir_marks = zeros(Int8, ld_cellsnum)
    d_bounds = spzeros(Int8, ld_cellsnum, 0)

    angles = Array{Array{Int64, 1}, 1}(undef, lld_cellsnum)

    for lld in 1:lld_cellsnum
        as = []
        for ld in ld_bounds[lld, :].nzind
            push!(as, (ld, angles_fn(lld, ld)))
        end
        sort!(as, lt=(a,b)->a[2]<b[2])
        as = map(a->a[1], as)
        angles[lld] = as
    end

    while (sigma = get_seed_cell(ld_cellsnum, count_marks)) > 0
        if verbose
            #print(Int(floor(50 * sum(count_marks) / ld_cellsnum)), "%\r")
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


function minimal_cycles(angles_fn::Function, verbose=false)
    return _minimal_cycles
end


EF = SparseArrays.sparse(ld_bounds)
