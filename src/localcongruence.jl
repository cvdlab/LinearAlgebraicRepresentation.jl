function chaincongruence(W, Delta_0, Delta_1)
@show W
@show Delta_0
@show Delta_1
    W = convert(Matrix, W')
    V, vclasses = vcongruence(W)
    EV, eclasses = cellcongruence(Delta_0, vclasses, dim=1)
    copEW = Lar.lar2cop(convert(Vector{Vector{Int64}}, EV))

    FE, fclasses = cellcongruence(Delta_1, eclasses, dim=2)
    copFE = Lar.lar2cop(convert(Vector{Vector{Int64}}, FE))
    
    copFV = convert(Lar.ChainOp, copFE * copEW .÷ 2)
    FV = convert(Lar.Cells, Lar.cop2lar(copFV))
    EV = convert(Lar.Cells, EV)
    copFE = Lar.build_copFE(FV::Lar.Cells, EV::Lar.Cells)
    copEV = Lar.coboundary_0(EV::Lar.Cells)
    
    V = convert(Matrix, V')
    return V, copEV, copFE
end

function cellcongruence(Delta, inclasses; dim)
    cellarray = Lar.cop2lar(Delta)
    newcell = Vector(undef, size(Delta,2))
    [ newcell[e] = k for (k, class) in enumerate(inclasses) for e in class ]
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

function vcongruence(V::Matrix; epsilon=1e-6)
   vclasses, visited = [], []
   kdtree = NearestNeighbors.KDTree(V);
   for vidx = 1 : size(V, 2) 
        if !(vidx in visited)
          nearvs = NearestNeighbors.inrange(kdtree, V[:,vidx], epsilon)
          push!(vclasses, nearvs)
          append!(visited, nearvs) 
        end
    end
   #W = hcat([sum(V[:,class], dims=2)/length(class) for class in vclasses]...)
@show vclasses
   sort!(map(sort!,vclasses))
@show vclasses
   W = hcat([V[:,class[1]] for class in vclasses]...)
   return W, vclasses
end
