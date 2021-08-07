"""
	vcongruence(
        V::Lar.Points; epsilon=1e-6
    )::Tuple{Lar.Points, Array{Array{Int,1},1}}

Evaluates the Vertex Congruence for 3D-points ``V ∈ ℳ(3,n)``.

The function determines the points of ``V`` closer than ``epsilon`` and
builds a new Vertex Set made of the representative of each point cluster.

The method returns:
 - the new Vertex Set
 - a map that, for every new vertex, identifies the old vertices it is made of
"""
function vcongruence(V::Lar.Points; epsilon=1e-6)
    vclasses, visited = [], []
    V = convert(Array{Float64,2}, V)
    kdtree = NearestNeighbors.KDTree(V)
    for vidx = 1 : size(V, 2)  if !(vidx in visited)
        nearvs = NearestNeighbors.inrange(kdtree, V[:,vidx], epsilon)
        push!(vclasses, nearvs)
        append!(visited, nearvs)
    end  end
    W = hcat([sum(V[:,class], dims=2)/length(class) for class in vclasses]...)
    return W, vclasses
end

"""
    cellCongruence(
        Delta::Lar.ChainOp,
        inclasses::Array{Array{Int,1},1};
        dim::Int = 0
    )::Tuple{Lar.Cells, Array{Array{Int,1},1}}

Evaluates the Cell Congruence for Cochain ``Delta`` with classes ``inclasses``.

The function determines the new Cochain (as ``Lar.Cells``) built from ``cop``
where the lower order cells are merged according to ``inclasses`` map.

The parameter ``dim``, if specified, represent the order of the cells
(that is the least number of lower order cells a current order cell is made of).

The method returns:
- the new Cochain Operator as ``Lar.Cells``
- a map that, for every new cell, identifies the old cells it is made of
"""

function cellcongruence(Delta, inclasses; dim = 0)
    cellarray = Lar.cop2lar(Delta)
    newcell = Vector(undef, size(Delta,2))
    for (k, class) in enumerate(inclasses)  for e in class
            newcell[e] = k
    end end
    cells = [sort(map(x -> newcell[x], face)) for face in cellarray]
    okcells = [cell for cell in cells if length(Set(cell)) > dim] # non-empty cells
    classes = DefaultOrderedDict{Vector, Vector}([])
    for (k,face) in enumerate(okcells)
        classes[face] == [] ? classes[face] = [k] : append!(classes[face], [k])
    end
    cells = collect(keys(classes))
    outclasses = collect(values(classes))
    return cells, outclasses
end

"""
    chaincongruence(
        W::Lar.Points,
        Delta_0::Lar.ChainOp,
        Delta_1::Lar.ChainOp
        epsilon = 1e-6
    )::Tuple{Lar.Points, Lar.Cells, Lar.Cells}

Cell Congruence Enabling Algorithm for 3D Complexes.

The algorithm performs the identify vertices within ``epsilon`` and coherently
update the Chain Operators.
"""
function chaincongruence(W, Delta_0, Delta_1; epsilon = 1e-6)
    V, vclasses = vcongruence(W, epsilon = epsilon)
    EV, eclasses = cellcongruence(Delta_0, vclasses, dim=1)
    FE, _ = cellcongruence(Delta_1, eclasses, dim=2)
    return V, EV, FE
end


"""
	signedCellCongruence(
		cop::Lar.ChainOp,
		lo_cls::Array{Array{Int,1},1},
		lo_sign::Array{Array{Int8,1},1};
		imp = false,
		d = 0
	)
Evaluates the Cell Congruence for a Cochain ``cop`` with classes ``lo_cls``.
The function determines the new Cochain Operator built from ``cop`` where
the lower order cells are merged according to ``lo_cls`` map.
``lo_sign`` specifies whether a cell must be considered in reverse order.
If optional paramenter ``imp`` is set to ``true`` then FP imprecisions
are taken into account in the sense that lower order cells may have collided.
The parameter ``d`` represent then the order of the cell (that also is the
least number of lower order cells a current order cell is made of).
The method returns:
 - the new Cochain Operator
 - a map that, for every new cell, identifies the old cells it is made of
 - a map that, for every new cell, specify if old cells have changed ordering.
"""
function signedCellCongruence(cop, lo_cls, lo_sign; imp = false, dim = 0)::Tuple{
        Lar.ChainOp,
		Array{Array{Int,1},1},
		Array{Array{Int8,1},1}
    }
#function cellCongruence(
#		cop::Lar.ChainOp,
#		lo_cls::Array{Array{Int,1},1},
#		lo_sign::Array{Array{Int8,1},1};
#		imp = false,
#		dim = 0
#	)::Tuple{
#		Lar.ChainOp,
#		Array{Array{Int,1},1},
#		Array{Array{Int8,1},1}
#	}

	# Collide columns
	copCols = []
	for i = 1 : length(lo_cls)
		col = sum([
			cop[:, lo_cls[i][j]] .* lo_sign[i][j]
			for j = 1 : length(lo_cls[i])
		])
#		if imp
#			#TODO remove zeros
#			if length(col.nzind) > dim  push!(copCols, col)  end
#		end
		push!(copCols, col)
	end

	# Retrieve Matrix Form and extract rows
	cop = hcat(copCols...)
	rows = [cop[row, :] for row = 1 : cop.m]

	# Adjustinfg signs such that the first value always is `-1`
	sign = ones(Int8, cop.m)
	for row = 1 : cop.m  if rows[row].nzval[1] > 0
		rows[row] = -rows[row]
		sign[row] = -1
	end  end

	# Sort rows in order to collide them quickly
	rows_ord = sortperm(rows)
	nrows = unique(rows[rows_ord])
	ho_cls = [Array{Int,1}() for i = 1 : length(nrows)]
	nidx = 1

	# Collide cells with same structure
	for cidx = 1 : cop.m
		if rows[rows_ord[cidx]] != nrows[nidx]  nidx = nidx + 1;  end
		push!(ho_cls[nidx], rows_ord[cidx])
	end

	# Reshaping Sign and Cochain
	ho_sign = [[sign[el] for el in cell] for cell in ho_cls]
#	cop = convert(Lar.ChainOp, hcat(nrows...)')

	return hcat(nrows...)', ho_cls, ho_sign
end