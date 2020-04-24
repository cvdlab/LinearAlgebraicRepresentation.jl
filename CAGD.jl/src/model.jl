import Base: +, length, size, ==, isempty

Lar = LinearAlgebraicRepresentation

#-------------------------------------------------------------------------------
#   MODEL STRUCT DEFINITION
#-------------------------------------------------------------------------------

mutable struct Model
    G::Lar.Points
    T::Array{Lar.ChainOp, 1}

    function Model(V::Lar.Points, T::Array{Lar.ChainOp, 1})
        dim, npts = size(V)
        dim > 0 ||
            throw(ArgumentError("At least one point is needed."))
        length(T) == dim ||
            throw(ArgumentError("Topology is not coherent with Geometry."))

        isempty(T[1]) || size(T[1], 2) == npts ||
            throw(ArgumentError("Topology not coherent with Geometry."))
        for i = 2 : dim
            isempty(T[i-1]) || isempty(T[i]) ||
                size(T[i], 2) == size(T[i-1],1) ||
                throw(ArgumentError("Topology not coherent with Topology."))
        end

        new(V, T)
    end

    function Model(V::Lar.Points)
        T = convert(Array{Lar.ChainOp,1},
            [SparseArrays.spzeros(Int8, 0, 0) for i = 1 : size(V, 1)]
        )
        T[1] = convert(Lar.ChainOp, SparseArrays.spzeros(Int8, 0, size(V,2)))
        CAGD.Model(V, T)
    end

    function Model()
        nothing
    end
end

#-------------------------------------------------------------------------------
#   BASIC PROPERTIES
#-------------------------------------------------------------------------------

length(m::CAGD.Model)               = size(m.G, 1)
size(m::CAGD.Model)		           = size(m.G)
size(m::CAGD.Model, i::Int)         = i == 0 ? size(m.G   ) : size(m.T[i]   )
size(m::CAGD.Model, i::Int, j::Int) = i == 0 ? size(m.G, j) : size(m.T[i], j)
==(m1::CAGD.Model, m2::CAGD.Model)   = m1.G == m2.G && m1.T == m2.T
isempty(m::CAGD.Model, d::Int)      = isempty(m.T[d])

#-------------------------------------------------------------------------------
#   OTHER BASE REIMPLEMENTATIONS
#-------------------------------------------------------------------------------

Base.copy(m::CAGD.Model)     = CAGD.Model(m.G, m.T)
Base.deepcopy(m::CAGD.Model) = CAGD.Model(Base.deepcopy(m.G), Base.deepcopy(m.T))

#-------------------------------------------------------------------------------
#   GEOMETRY MANIPULATION
#-------------------------------------------------------------------------------

function addModelVertices!(m::CAGD.Model, V::Lar.Points)::Nothing
    length(m) == size(V, 1) ||
        throw(ArgumentError("Point dimension mismatch."))
    m.G = [m.G V]

    m.T[1] = SparseArrays.sparse(
        findnz(m.T[1])...,
        size(m, 1, 1),
        size(m, 1, 2) + size(V, 2)
    )
    return
end

function addModelVertex!(m::CAGD.Model, v::Array{Float64, 1})::Nothing
    return addModelVertices!(m, v[:, :])
end

function deleteModelVertex!(m::CAGD.Model, v::Int)::Nothing
    deleteModelVertices!(m, [v])
end

function deleteModelVertices!(m::CAGD.Model, vs::Array{Int,1})::Nothing
    tokeepV = setdiff(collect(1 : size(m, 0, 2)), vs)
    m.G     = m.G[:, tokeepV]

    if !isempty(m.T[1])
        # From EV must be deleted all the columns linked to a deleted vertex
        #  and all rows related to a dangling edge
        todelE = abs.(m.T[1][:, vs])
        todelE = [l for l = 1 : todelE.m if sum(todelE[l, :]) != 0]
        m.T[1] = m.T[1][:, tokeepV]
        if !isempty(todelE) CAGD.deleteModelCells!(m, 1, todelE) end
    end

    return
end

#-------------------------------------------------------------------------------
#   TOPOLOGY MANIPULATION
#-------------------------------------------------------------------------------

function addModelCells!(m::CAGD.Model, deg::Int, cs::Lar.ChainOp)::Nothing
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    (cs.n == size(m, deg - 1, 1) && deg > 1) ||
        (cs.n == size(m, deg - 1, 2) && deg == 1) ||
        throw(ArgumentError("Incoherent Chain Operator size"))

    # Adding input cells to the topological structure
    m.T[deg] = [m.T[deg]; cs]

    # Adding slot for the new cells to the higher order cells
    if deg < length(m)
        m.T[deg + 1] = [
            m.T[deg + 1] SparseArrays.spzeros(Int8, size(m, deg + 1, 1), cs.m)
        ]
    end

    return
end

function addModelCell!(m::CAGD.Model, deg::Int, c::Lar.Cell)::Nothing
    CAGD.addModelCells!(m, deg, convert(Lar.ChainOp, c))
end

function addModelCells!(
        m::CAGD.Model, deg::Int, cs::Lar.Cells; signed=false
    )::Nothing

    if signed
        if deg == 1
            scs = convert(Lar.ChainOp, Lar.coboundary_0(cs))
        elseif deg == 2
            EV = Lar.cop2lar(m.T[1]);
            FV = [collect(Set(vcat([EV[e] for e in face]...))) for face in cs]
            scs = convert(Lar.ChainOp, Lar.coboundary_1(m.G, FV, EV))
        else
            println("No Methods for $deg-signed-coboundary.")
            println("Using unsigned coboundary instead.")
            signed = false
        end
    end

    if !signed
        I = Array{Int,1}()
        J = Array{Int,1}()
        K = Array{Int8,1}()
        for i = 1 : length(cs)
            for j = 1 : length(cs[i])
                push!(I, i)
                push!(J, cs[i][j])
                push!(K, 1)
            end
        end

        scs = SparseArrays.sparse(I, J, K, length(cs), size(m, deg, 2))
    end

    return CAGD.addModelCells!(m, deg, scs)
end

function addModelCell!(m::CAGD.Model, deg::Int, c::Array{Int64,1})::Nothing
    CAGD.addModelCells!(m, deg, [c])
end


function deleteModelCells!(m::CAGD.Model, deg::Int, cs::Array{Int, 1})::Nothing
    !isempty(cs) || return
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    !isempty(m.T[deg]) ||
        throw(ArgumentError("There are no cells of degree $deg"))
    max(cs...) ≤ size(m, deg, 1) ||
        throw(ArgumentError("There are not $(max(cs...)) cells of degree $deg"))

    # Removing `cs` rows from `m.T[deg]`
    tokeep     = setdiff(collect(1 : size(m, deg, 1)), cs)
    m.T[deg]   = m.T[deg][tokeep, :]

    if deg == length(m) return end

    # Removing `cs` cols from `m.T[deg+1]` by checking if some `deg+1` cell
    #  has to be removed too.
    # An higher ord cell must be removed if it contains a removed lower ord cell
    todelHo    = m.T[deg + 1][:, cs]
    m.T[deg+1] = m.T[deg + 1][:, tokeep]
    if isempty(m.T[deg+1]) return end
    todelHo    = [l for l = 1 : todel.m if sum(todelE[l, :]) != 0]
    if !isempty(todelHo) CAGD.deleteModelCells!(m, deg+1, todelHo) end
    return
end

function deleteModelCell!(m::CAGD.Model, deg::Int, c::Int)::Nothing
    CAGD.deleteModelCells!(m, deg, [c])
end

function getModelLoCell(m::CAGD.Model, deg::Int, c::Int)::Array{Int, 1}
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    return m.T[deg][c, :].nzind
end

function getModelCellVertices(m::CAGD.Model, deg::Int, c::Int, ret_idx=false)
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    set = [c]
    for d = deg : -1 : 1
        set = ∪([m.T[d][el, :].nzind for el in set]...)
    end

    if ret_idx == true
        return (map(i -> m.G[:, set[i]], 1:length(set)), set)
    end
    return map(i -> m.G[:, set[i]], 1:length(set))
end

function getModelCellGeometry(m::CAGD.Model, deg::Int, c::Int, ret_idx=false)
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    set = [c]
    for d = deg : -1 : 1
        set = ∪([m.T[d][el, :].nzind for el in set]...)
    end

    if ret_idx == true
        return (m.G[:, set], set)
    end
    return m.G[:, set]
end

#-------------------------------------------------------------------------------
#   MODEL MANIPULATION
#-------------------------------------------------------------------------------

"""
    modelPurge!(m::CAGD.Model, [depth::Int = 0])::Nothing

Purges the model Geometry from unused elements.

Purges each Geometry element that is not related to a topological Cell.
If the optional argument `depth` is specified than it also purges all the
topological cells of degree lower or equal to depth that are not related
to an higher order topological cell.
"""
function modelPurge!(m::CAGD.Model, depth::Int = 0)::Nothing
    depth < length(m) ||
        throw(ArgumentError("The model do not have degree $(depth+1)"))

    for d = depth : -1 : 1
        todel = [i for i = 1 : size(m, d+1, 2) if sum(abs.(m.T[d+1][:,i])) == 0]
        if !isempty(todel) CAGD.deleteModelCells!(m, d, todel) end
    end
    todel = [i for i = 1 : size(m, 1, 2) if sum(abs.(m.T[1][:, i])) == 0]
    if !isempty(todel) CAGD.deleteModelVertices!(m, todel) end
    return
end


#-------------------------------------------------------------------------------
#   MODEL INTERACTION
#-------------------------------------------------------------------------------

function uniteModels!(model::CAGD.Model, m2::CAGD.Model)::Nothing

    length(model) == length(m2) ||
        throw(ArgumentError("ERROR: Inconsistent models dimension!"))

    CAGD.addModelVertices!(model, m2.G)
    for d = 1 : length(model)
        shift = size(model, d, 2) - size(m2, d, 2)
        I, J, K = SparseArrays.findnz(m2.T[d])
        J = J .+ shift
        CAGD.addModelCells!(model, d, SparseArrays.sparse(
            I, J, K, size(m2, d, 1), size(m2, d, 2) + shift
        ))
    end

    return
end

function uniteModels(m1::CAGD.Model, m2::CAGD.Model)::CAGD.Model
    model = CAGD.deepcopy(m1)
    CAGD.uniteModels!(model, m2)
    return model
end

function uniteMultipleModels(models::Array{CAGD.Model,1})::CAGD.Model

    model = deepcopy(models[1])

    for i = 2 : length(models)
        CAGD.uniteModels!(model, models[i])
    end

    return model
end

function mergeMultipleModels(models::Array{CAGD.Model,1}; err=1e-6)::CAGD.Model
    return CAGD.mergeModelVertices(uniteMultipleModels(models), err=err)
end

function mergeModelVertices(model::CAGD.Model; err=1e-6, signed_merge = false)
    V, cls = CAGD.vcongruence(model.G, epsilon = err)
    if signed_merge  lo_sign = [ones(Int8, length(cl)) for cl in cls]  end
    congModel = CAGD.Model(V)
    for d = 1 : length(model)
        if isempty(model.T[d])  break  end
        if signed_merge
            cop, cls, lo_sign = CAGD.signedCellCongruence(model.T[d], cls, lo_sign, dim=d)
            if  d > 1 cop = -cop  end
        else
            cop, cls = CAGD.cellcongruence(model.T[d], cls, dim=d)
            cop = convert(Lar.Cells, cop)
        end
        CAGD.addModelCells!(congModel, d, cop)
    end

    return congModel
end

#-------------------------------------------------------------------------------
#   MODEL EXPORT
#-------------------------------------------------------------------------------

function jlexportModel(
        model::CAGD.Model,
        filename::String,
        modelname::String = "model"
    )::Nothing

    open(filename, "w") do f
        # Write imports
        write(f, "# using LinearAlgebraicRepresentation\n")
        write(f, "using SparseArrays\n")
        write(f, "# using CAGD\n\n")

        # Write Geometry
        write(f, "#=== GEOMETRY ===#\n")
        write(f, "$modelname = CAGD.Model([")
        for d = 1 : length(model)
            write(f, "\n\t")
            for p in model.G[d, :]  write(f, "$p ")  end
        end
        write(f, "\n]);\n\n")

        # Write Topology
        for d = 1 : length(model)
            write(f, "#=== TOPOLOGY : $d-CELLS ===#\n")
            I, J, K = findnz(model.T[d])
            write(f, "I = Array{Int64}([")
            for i in I  write(f, "$i, ")  end
            write(f, "]);\nJ = Array{Int64}([")
            for j in J  write(f, "$j, ")  end
            write(f, "]);\nK = Array{Int8}([")
            for k in K  write(f, "$k, ")  end
            write(f, "]);\n\n")
            write(f, "CAGD.addModelCells!($modelname, $d, SparseArrays.sparse(")
            write(f, "I, J, K, $(size(model, d, 1)), $(size(model, d, 2))")
            write(f, "))\n\n\n")
        end
    end

    return
end
