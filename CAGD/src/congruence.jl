using LinearAlgebraicRepresentation
using NearestNeighbors

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
    cells = [map(x -> newcell[x], face) for face in cellarray]
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
