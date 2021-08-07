import Base: +, length, size, ==, isempty, show

Lar = LinearAlgebraicRepresentation

#-------------------------------------------------------------------------------
#   MODEL STRUCT DEFINITION
#-------------------------------------------------------------------------------

mutable struct Model
    G::Lar.Points
    T::Array{Lar.ChainOp, 1}

    """
        Model(V::Lar.Points, T::Array{Lar.ChainOp, 1})
    
    Generic constructor for CAGD.Model.
    Coherency checks are performed between Topology and Geometry
    """
    function Model(V::Lar.Points, T::Array{Lar.ChainOp, 1})
        dim, npts = size(V)
        dim > 0 ||
            throw(ArgumentError("At least one point is needed."))
        length(T) == dim ||
            throw(ArgumentError("Topology is not coherent with Geometry."))

        size(T[1], 2) == npts ||
            throw(ArgumentError("Topology not coherent with Geometry."))
        for i = 2 : dim
            isempty(T[i-1]) || size(T[i], 2) == size(T[i-1],1) ||
                throw(ArgumentError("Topology not coherent with Topology."))
        end

        new(V, T)
    end

    function Model(V::Lar.Points, T::Array{Lar.Cells, 1})
        m = Model(V)
        addModelCells!(m, 1, T[1], signed=true)
        cFE = convert(Lar.ChainOp, Lar.coboundary_1(m.G, T[2], T[1]))
        addModelCells!(m, 2, cFE)
        return m
    end

    """
        Model(V::Lar.Points, EV::Lar.Cells)
    
    Builds a signed CAGD.Model with vertices and edges only.
    """
    function Model(V::Lar.Points, EV::Lar.Cells)

        I = Array{Int,1}()
        J = Array{Int,1}()
        K = Array{Int8,1}()
        for i = 1 : length(EV)
            let sign = -1;  for j = 1 : length(EV[i])
                push!(I, i)
                push!(J, EV[i][j])
                push!(K, sign)
                sign = 1
            end  end
        end

        T = [SparseArrays.sparse(I, J, K, length(EV), size(V, 2))]
        if size(V, 1) > 1  push!(T, SparseArrays.spzeros(Int8, 0, length(EV)))  end
        Ts = SparseArrays.spzeros(Int8, 0, 0)
        for i = 3 : size(V, 1)  push!(T, Ts)  end
        CAGD.Model(V, T)
    end

    """
        Model(V::Lar.Points)
    
    Constructor for CAGD.Model with geometry only.
    Topology is set to void by default.
    """
    function Model(V::Lar.Points)
        T = convert(Array{Lar.ChainOp,1},
            [SparseArrays.spzeros(Int8, 0, 0) for i = 1 : size(V, 1)]
        )
        T[1] = convert(Lar.ChainOp, SparseArrays.spzeros(Int8, 0, size(V,2)))
        CAGD.Model(V, T)
    end

    """
        Model()
    
    Void constructor for CAGD.Model.
    Returns `nothing`
    """
    function Model()
        nothing
    end
end

#-------------------------------------------------------------------------------
#   BASIC PROPERTIES
#-------------------------------------------------------------------------------

length(m::CAGD.Model)               = size(m.G, 1)
size(m::CAGD.Model)		            = size(m.G)
size(m::CAGD.Model, i::Int)         = i == 0 ? size(m.G   ) : size(m.T[i]   )
size(m::CAGD.Model, i::Int, j::Int) = i == 0 ? size(m.G, j) : size(m.T[i], j)
==(m1::CAGD.Model, m2::CAGD.Model)  = m1.G == m2.G && m1.T == m2.T
isempty(m::CAGD.Model, d::Int)      = isempty(m.T[d])
function show(io::IO, m::CAGD.Model)
    println(io, "$(length(m))D model with:")
    println(io, " - $(size(m, 1, 2)) points")
    for d = 1 : length(m)  println(io, " - $(size(m, d, 1)) $d-cells")  end
    display(m.G)
    for d = 1 : length(m)  display(Matrix(m.T[d]))  end
end

#-------------------------------------------------------------------------------
#   OTHER BASE REIMPLEMENTATIONS
#-------------------------------------------------------------------------------

Base.copy(m::CAGD.Model)     = CAGD.Model(m.G, m.T)
Base.deepcopy(m::CAGD.Model) = CAGD.Model(Base.deepcopy(m.G), Base.deepcopy(m.T))

#-------------------------------------------------------------------------------
#   GEOMETRY MANIPULATION
#-------------------------------------------------------------------------------

"""
    addModelVertices!(m::CAGD.Model, V::Lar.Points)::Nothing

Adds `V` vertices to model `m` in tail of `m.G`.
Coherently updates `m.T[1]`.
Points dimension must match.

---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0;0.0 1.0]);

julia> CAGD.addModelVertices!(m, [1.0 1.0; 0.0 1.0])

julia> m
2D model with:
 - 4 points
 - 0 1-cells
 - 0 2-cells
2×4 Array{Float64,2}:
 0.0  0.0  1.0  1.0
 0.0  1.0  0.0  1.0
0×4 Array{Int8,2}
0×0 Array{Int8,2}
```
"""
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

"""
    addModelVertex!(m::CAGD.Model, v::Array{Float64,1})::Nothing

Adds `v` vertex to model `m` in tail of `m.G`.
Coherently updates `m.T[1]`.
Point dimension must match.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0;0.0 1.0]);

julia> CAGD.addModelVertex!(m, [1.0; 1.0])

julia> m
2D model with:
 - 3 points
 - 0 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  1.0
0×3 Array{Int8,2}
0×0 Array{Int8,2}
```
"""
function addModelVertex!(m::CAGD.Model, v::Array{Float64, 1})::Nothing
    return addModelVertices!(m, v[:, :])
end

"""
    deleteModelVertices!(m::CAGD.Model, vs::Array{Int,1})::Nothing

Delete `vs` indexed vertices from model `m`.
Coherently deletes cells containing such points.
---
# Examples
```jldoctest
julia>  m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
           SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
           SparseArrays.spzeros(Int8, 0, 3)
       ]);

julia> CAGD.deleteModelVertices!(m, [1; 2])

julia> m
2D model with:
 - 1 points
 - 0 1-cells
 - 0 2-cells
2×1 Array{Float64,2}:
 1.0
 0.0
0×1 Array{Int8,2}
0×3 Array{Int8,2}
```
"""
function deleteModelVertices!(m::CAGD.Model, vs::Array{Int,1})::Nothing
    min(vs...) > 0 || throw(ArgumentError("Vertex indices are positive integers"))
    max(vs...) ≤ size(m, 0, 2) || throw(ArgumentError("There are not $(max(vs...)) vertices"))

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

"""
    deleteModelVertex!(m::CAGD.Model, v::Int)::Nothing

Delete `v` indexed vertex from model `m`.
Coherently deletes cells containing such point.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.spzeros(Int8, 0, 3)
]);

julia> CAGD.deleteModelVertex!(m, 1)

julia> m
2D model with:
 - 2 points
 - 1 1-cells
 - 0 2-cells
2×2 Array{Float64,2}:
 0.0  1.0
 1.0  0.0
1×2 Array{Int8,2}:
 -1  1
0×3 Array{Int8,2}
```
"""
function deleteModelVertex!(m::CAGD.Model, v::Int)::Nothing
    deleteModelVertices!(m, [v])
end

#-------------------------------------------------------------------------------
#   TOPOLOGY MANIPULATION
#-------------------------------------------------------------------------------

"""
    addModelCells!(m::CAGD.Model, deg::Int, cs::Lar.ChainOp)::Nothing

Adds to model `m` the `deg` cells specified by `cs` given as sparse matrix.
Coherently updates `deg+1` cells.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0]);

julia> CAGD.addModelCells!(m, 1, SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]))

julia> m
2D model with:
 - 3 points
 - 3 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
3×3 Array{Int8,2}:
 -1   1  0
 -1   0  1
  0  -1  1
0×3 Array{Int8,2}
```
"""
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

"""
    addModelCell!(m::CAGD.Model, deg::Int, c::Lar.Cell)::Nothing

Adds to model `m` the `deg` cell specified by `c` as sparse vector.
Coherently updates `deg+1` cells.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0]);

julia> CAGD.addModelCell!(m, 1, SparseArrays.sparse(Int8[-1; 1; 0]))

julia> m
2D model with:
 - 3 points
 - 1 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
1×3 Array{Int8,2}:
 -1  1  0
0×1 Array{Int8,2}
```
"""
function addModelCell!(m::CAGD.Model, deg::Int, c::Lar.Cell)::Nothing
    CAGD.addModelCells!(m, deg, convert(Lar.ChainOp, c'))
end

"""
    addModelCells!(m::CAGD.Model, deg::Int, cs::Lar.Cells; signed=false)::Nothing

Adds to model `m` the `deg` cells specified by `cs` as array of arrays.
Coherently updates `deg+1` cells.
If `signed` is set to `true`, cells are evaluated with sign.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0]);

julia> CAGD.addModelCells!(m, 1, [[1, 2], [1, 3], [2, 3]])

julia> CAGD.addModelCells!(m, 1, [[1, 2], [1, 3], [2, 3]], signed = true)

julia> m
2D model with:
 - 3 points
 - 6 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
6×3 Array{Int8,2}:
  1   1  0
  1   0  1
  0   1  1
 -1   1  0
 -1   0  1
  0  -1  1
0×6 Array{Int8,2}
```
"""
function addModelCells!(
        m::CAGD.Model, deg::Int, cs::Lar.Cells; signed=false
    )::Nothing

    if signed
        if deg == 1
            I = Array{Int,1}()
            J = Array{Int,1}()
            K = Array{Int8,1}()
            for i = 1 : length(cs)
                let sign = -1;  for j = 1 : length(cs[i])
                    push!(I, i)
                    push!(J, cs[i][j])
                    push!(K, sign)
                    sign = 1
                end  end
            end

            scs = SparseArrays.sparse(I, J, K, length(cs), size(m, deg, 2))
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

"""
    addModelCell!(m::CAGD.Model, deg::Int, c::Array{Int64,1}; signed=false)::Nothing

Adds to model `m` the `deg` cell specified by `c` as index of `deg-1` cells.
Coherently updates `deg+1` cells.
If `signed` is set to `true`, the cell is computed with sign.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0]);

julia> CAGD.addModelCell!(m, 1, [1, 2])

julia> CAGD.addModelCell!(m, 1, [1, 3], signed = true)

julia> m
2D model with:
 - 3 points
 - 2 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
2×3 Array{Int8,2}:
  1  1  0
 -1  0  1
0×2 Array{Int8,2}
```
"""
function addModelCell!(m::CAGD.Model, deg::Int, c::Array{Int64,1}; signed=false)::Nothing
    CAGD.addModelCells!(m, deg, [c], signed = signed)
end

"""
    deleteModelCells!(m::CAGD.Model, deg::Int, cs::Array{Int, 1})::Nothing

Delete `cs` indexed `deg`-cells from model `m`.
Coherently deletes higher order cells containing such cells.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> CAGD.deleteModelCells!(m, 1, [1; 2])

julia> m
2D model with:
 - 3 points
 - 1 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
1×3 Array{Int8,2}:
 0  -1  1
0×1 Array{Int8,2}
```
"""
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

    if  deg == length(m) return  end

    # Removing `cs` cols from `m.T[deg+1]` by checking if some `deg+1` cell
    #  has to be removed too.
    # An higher ord cell must be removed if it contains a removed lower ord cell
    todelHo    = m.T[deg + 1][:, cs]
    m.T[deg+1] = m.T[deg + 1][:, tokeep]
    if isempty(m.T[deg+1]) return end
    todelHo, _, _ = SparseArrays.findnz(todelHo)
    if !isempty(todelHo) CAGD.deleteModelCells!(m, deg+1, unique(todelHo)) end
    return
end

"""
    deleteModelCells!(m::CAGD.Model, deg::Int, c::Int)::Nothing

Delete `c` indexed `deg`-cell from model `m`.
Coherently deletes higher order cells containing such cell.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> CAGD.deleteModelCell!(m, 1, 1)

julia> m
2D model with:
 - 3 points
 - 2 1-cells
 - 0 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
2×3 Array{Int8,2}:
 -1   0  1
  0  -1  1
0×2 Array{Int8,2}
```
"""
function deleteModelCell!(m::CAGD.Model, deg::Int, c::Int)::Nothing
    CAGD.deleteModelCells!(m, deg, [c])
end

"""
    getModelLoCell(m::CAGD.Model, deg::Int, c::Int)::Array{Int, 1}
    getModelLoCell(m::CAGD.Model, deg::Int, cs::Array{Int,1})::Array{Int, 1}

Determines the indices of `deg-1` cells forming the boundary of `c`/`cs` indexed `deg` cell(s) of `m`
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0], [
    SparseArrays.sparse(Int8[-1 1 0 0; -1 0 1 0; 0 -1 1 0; 0 -1 0 1; 0 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1 0 0; 0 0 1 -1 1])
])
2D model with:
 - 4 points
 - 5 1-cells
 - 2 2-cells
2×4 Array{Float64,2}:
 0.0  0.0  1.0  1.0
 0.0  1.0  0.0  1.0
5×4 Array{Int8,2}:
 -1   1   0  0
 -1   0   1  0
  0  -1   1  0
  0  -1   0  1
  0   0  -1  1
2×5 Array{Int8,2}:
 1  -1  1   0  0
 0   0  1  -1  1

julia> CAGD.getModelLoCell(m, 2, 2)
3-element Array{Int64,1}:
 3
 4
 5

julia> CAGD.getModelLoCell(m, 1, [3, 4, 5])
3-element Array{Int64,1}:
 2
 3
 4
```
"""
function getModelLoCell(m::CAGD.Model, deg::Int, c::Int)::Array{Int, 1}
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    c > 0 || throw(ArgumentError("Cell indices are positive integers"))
    c ≤ size(m, deg, 1) || throw(ArgumentError("There are not $c $deg-cells"))
    return m.T[deg][c, :].nzind
end

function getModelLoCell(m::CAGD.Model, deg::Int, cs::Array{Int, 1})::Array{Int, 1}
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    min(cs...) > 0 || throw(ArgumentError("Cell indices are positive integers"))
    max(cs...) ≤ size(m, deg, 1) || throw(ArgumentError("There are not $(max(cs...)) $deg-cells"))
    return ∪([m.T[deg][c, :].nzind for c in cs]...)
end

"""
    getModelCellVertices(m::CAGD.Model, deg::Int, c::Int, ret_idx=false)
    getModelCellVertices(m::CAGD.Model, deg::Int, cs::Array{Int,1}, ret_idx=false)
        ::Union{Array{Array{Float64,1},1}, Tuple{Array{Array{Float64,1},1}, Array{Int,1}}}

Determines the list of vertices of `c`/`cs` indexed `deg` cell(s) of `m`.
If `ret_idx` is set true, it also returns the vertices list of index.

See also [`CAGD.getModelCellGeometry`](@ref) to get vertices as `Lar.Points`.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0], [
    SparseArrays.sparse(Int8[-1 1 0 0; -1 0 1 0; 0 -1 1 0; 0 -1 0 1; 0 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1 0 0; 0 0 1 -1 1])
])
2D model with:
 - 4 points
 - 5 1-cells
 - 2 2-cells
2×4 Array{Float64,2}:
 0.0  0.0  1.0  1.0
 0.0  1.0  0.0  1.0
5×4 Array{Int8,2}:
 -1   1   0  0
 -1   0   1  0
  0  -1   1  0
  0  -1   0  1
  0   0  -1  1
2×5 Array{Int8,2}:
 1  -1  1   0  0
 0   0  1  -1  1

julia> CAGD.getModelCellVertices(m, 2, 2)
3-element Array{Array{Float64,1},1}:
 [0.0, 1.0]
 [1.0, 0.0]
 [1.0, 1.0]

julia> CAGD.getModelCellVertices(m, 1, [3, 4, 5], true)
(Array{Float64,1}[[0.0, 1.0], [1.0, 0.0], [1.0, 1.0]], [2, 3, 4])
```
"""
function getModelCellVertices(
        m::CAGD.Model, deg::Int, cs::Array{Int,1}, ret_idx=false
    )::Union{Array{Array{Float64,1},1}, Tuple{Array{Array{Float64,1},1}, Array{Int,1}}}
    
    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    min(cs...) > 0 || throw(ArgumentError("Cell indices are positive integers"))
    max(cs...) ≤ size(m, deg, 1) || throw(ArgumentError("There are not $(max(cs...)) $deg-cells"))

    for d = deg : -1 : 1
        cs = ∪([m.T[d][el, :].nzind for el in cs]...)
    end

    if ret_idx == true
        return (map(i -> m.G[:, cs[i]], 1:length(cs)), cs)
    end
    return map(i -> m.G[:, cs[i]], 1:length(cs))
end

function getModelCellVertices(
        m::CAGD.Model, deg::Int, c::Int, ret_idx=false
    )::Union{Array{Array{Float64,1},1}, Tuple{Array{Array{Float64,1},1}, Array{Int,1}}}
    return getModelCellVertices(m, deg, [c], ret_idx)
end

"""
    getModelCellGeometry(m::CAGD.Model, deg::Int, c::Int, ret_idx=false)
    getModelCellGeometry(m::CAGD.Model, deg::Int, cs::Array{Int,1}, ret_idx=false)
        ::Union{Lar.Points, Tuple{Lar.Points, Array{Int,1}}}

Determines the geometry of `c`/`cs` indexed `deg` cell(s) of `m`.
If `ret_idx` is set true, it also returns the vertices list of index.

See also [`CAGD.getModelCellVertices`](@ref) to get vertices as a list.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0 1.0; 0.0 1.0 0.0 1.0], [
    SparseArrays.sparse(Int8[-1 1 0 0; -1 0 1 0; 0 -1 1 0; 0 -1 0 1; 0 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1 0 0; 0 0 1 -1 1])
])
2D model with:
 - 4 points
 - 5 1-cells
 - 2 2-cells
2×4 Array{Float64,2}:
 0.0  0.0  1.0  1.0
 0.0  1.0  0.0  1.0
5×4 Array{Int8,2}:
 -1   1   0  0
 -1   0   1  0
  0  -1   1  0
  0  -1   0  1
  0   0  -1  1
2×5 Array{Int8,2}:
 1  -1  1   0  0
 0   0  1  -1  1

julia> CAGD.getModelCellGeometry(m, 2, 2)
2×3 Array{Float64,2}:
 0.0  1.0  1.0
 1.0  0.0  1.0

julia> CAGD.getModelCellGeometry(m, 1, [3, 4, 5], true)
([0.0 1.0 1.0; 1.0 0.0 1.0], [2, 3, 4])
```
"""
function getModelCellGeometry(
        m::CAGD.Model, deg::Int, cs::Array{Int,1}, ret_idx=false
    )::Union{Lar.Points, Tuple{Lar.Points, Array{Int,1}}}

    deg > 0 || throw(ArgumentError("Degree must be a non negative value"))
    deg ≤ length(m) || throw(ArgumentError("The model do not have degree $deg"))
    min(cs...) > 0 || throw(ArgumentError("Cell indices are positive integers"))
    max(cs...) ≤ size(m, deg, 1) || throw(ArgumentError("There are not $(max(cs...)) $deg-cells"))
    for d = deg : -1 : 1
        cs = ∪([m.T[d][el, :].nzind for el in cs]...)
    end

    if ret_idx == true
        return (m.G[:, cs], cs)
    end
    return m.G[:, cs]
end

function getModelCellGeometry(
        m::CAGD.Model, deg::Int, c::Int, ret_idx=false
    )::Union{Lar.Points, Tuple{Lar.Points, Array{Int,1}}}

    return getModelCellGeometry(m, deg, [c], ret_idx)
end


#-------------------------------------------------------------------------------
#   MODEL MANIPULATION
#-------------------------------------------------------------------------------

"""
    modelPurge!(m::CAGD.Model[; depth::Int = 0])::Nothing

Purges the model Geometry from unused elements.

Purges each Geometry element that is not related to a topological Cell.
If the optional argument `depth` is specified than it also purges all the
topological cells of degree lower or equal to `depth` that are not related
to an higher order topological cell.
---
# Examples
```jldoctest
julia> m = CAGD.Model([0.0 0.0 1.0 0.0 2.0; 0.0 1.0 0.0 2.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0 0 0; -1 0 1 0 0; 0 -1 1 0 0; 0 -1 0 1 0]),
    SparseArrays.sparse(Int8[1 -1 1 0 ])
])
2D model with:
 - 5 points
 - 4 1-cells
 - 1 2-cells
2×5 Array{Float64,2}:
 0.0  0.0  1.0  0.0  2.0
 0.0  1.0  0.0  2.0  0.0
4×5 Array{Int8,2}:
 -1   1  0  0  0
 -1   0  1  0  0
  0  -1  1  0  0
  0  -1  0  1  0
1×4 Array{Int8,2}:
 1  -1  1  0

julia> CAGD.modelPurge!(m)

julia> m
2D model with:
 - 4 points
 - 4 1-cells
 - 1 2-cells
2×4 Array{Float64,2}:
 0.0  0.0  1.0  0.0
 0.0  1.0  0.0  2.0
4×4 Array{Int8,2}:
 -1   1  0  0
 -1   0  1  0
  0  -1  1  0
  0  -1  0  1
1×4 Array{Int8,2}:
 1  -1  1  0

julia> CAGD.modelPurge!(m, depth = 1)

julia>m
2D model with:
 - 3 points
 - 3 1-cells
 - 1 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
3×3 Array{Int8,2}:
 -1   1  0
 -1   0  1
  0  -1  1
1×3 Array{Int8,2}:
 1  -1  1
```
"""
function modelPurge!(m::CAGD.Model; depth::Int = 0)::Nothing
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

"""
    mergeModelVertices(m::CAGD.Model[; err=1e-6[, signed_merge = false]])::CAGD.Model

Evaluates Geometry congruence with atol `err` and coherently updates Topology via CCE algorithm of `m`.

The method regroups vertices closer than `err` in a single vertex via KD-tree.
Then it coherently updates topology cells, eventually mashing them up if needed.
A topological gift wrapping (see [`CAGD.tgw`](@ref)) might be needed after the application,
if cells merged loosing in resolution.

If `signed_merge` is set to true, than `m` is expected signed and it returns a signed model.

See also:
 - [`CAGD.vcongruence`](@ref)
 - [`CAGD.cellcongruence`](@ref)
 - [`CAGD.signedCellCongruence`](@ref)
---
# Examples
```jldoctest

```
"""
function mergeModelVertices(m::CAGD.Model; err=1e-6, signed_merge = false)
    V, cls = CAGD.vcongruence(m.G, epsilon = err)
    if signed_merge  lo_sign = [ones(Int8, length(cl)) for cl in cls]  end
    congModel = CAGD.Model(V)
    for d = 1 : length(m)
        if isempty(m.T[d])  break  end
        if signed_merge
            cop, cls, lo_sign = CAGD.signedCellCongruence(m.T[d], cls, lo_sign, dim=d)
            if  d > 1 cop = -cop  end
        else
            cop, cls = CAGD.cellcongruence(m.T[d], cls, dim=d)
            cop = convert(Lar.Cells, cop)
        end
        CAGD.addModelCells!(congModel, d, cop)
    end

    return congModel
end


#-------------------------------------------------------------------------------
#   MODEL INTERACTION
#-------------------------------------------------------------------------------

"""
    uniteModels!(model::CAGD.Model, m2::CAGD.Model)::Nothing

Fuses `m2` in `model` by adding coherently geometry and topology.
---
# Examples
```jldoctest
julia> m1 = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> m2 = CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> CAGD.uniteModels!(m1, m2)

julia> m1
2D model with:
 - 6 points
 - 6 1-cells
 - 2 2-cells
2×6 Array{Float64,2}:
 0.0  0.0  1.0  0.0   0.0  -1.0
 0.0  1.0  0.0  0.0  -1.0   0.0
6×6 Array{Int8,2}:
 -1   1  0   0   0  0
 -1   0  1   0   0  0
  0  -1  1   0   0  0
  0   0  0  -1   1  0
  0   0  0  -1   0  1
  0   0  0   0  -1  1
2×6 Array{Int8,2}:
 1  -1  1  0   0  0
 0   0  0  1  -1  1
```
"""
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

"""
    uniteModels(m1::CAGD.Model, m2::CAGD.Model)::CAGD.Model

Fuses `m1` and `m2` in a new model by adding coherently geometry and topology.
---
# Examples
```jldoctest
julia> m1 = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> m2 = CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> CAGD.uniteModels(m1, m2)
2D model with:
 - 6 points
 - 6 1-cells
 - 2 2-cells
2×6 Array{Float64,2}:
 0.0  0.0  1.0  0.0   0.0  -1.0
 0.0  1.0  0.0  0.0  -1.0   0.0
6×6 Array{Int8,2}:
 -1   1  0   0   0  0
 -1   0  1   0   0  0
  0  -1  1   0   0  0
  0   0  0  -1   1  0
  0   0  0  -1   0  1
  0   0  0   0  -1  1
2×6 Array{Int8,2}:
 1  -1  1  0   0  0
 0   0  0  1  -1  1
```
"""
function uniteModels(m1::CAGD.Model, m2::CAGD.Model)::CAGD.Model
    model = CAGD.deepcopy(m1)
    CAGD.uniteModels!(model, m2)
    return model
end

"""
    uniteMultipleModels(models::Array{CAGD.Model,1})::CAGD.Model

Fuses all the `models` in a new model by adding coherently Geometry and Topology
---
# Examples
```jldoctest
julia> models = [
    CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ]);

    CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ]);
    
    CAGD.Model([0.0 -1.0 -1.0; 1.0 0.0 1.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ])
];

julia> CAGD.uniteMultipleModels(models)
2D model with:
 - 9 points
 - 9 1-cells
 - 3 2-cells
2×9 Array{Float64,2}:
 0.0  0.0  1.0  0.0   0.0  -1.0  0.0  -1.0  -1.0
 0.0  1.0  0.0  0.0  -1.0   0.0  1.0   0.0   1.0
9×9 Array{Int8,2}:
 -1   1  0   0   0  0   0   0  0
 -1   0  1   0   0  0   0   0  0
  0  -1  1   0   0  0   0   0  0
  0   0  0  -1   1  0   0   0  0
  0   0  0  -1   0  1   0   0  0
  0   0  0   0  -1  1   0   0  0
  0   0  0   0   0  0  -1   1  0
  0   0  0   0   0  0  -1   0  1
  0   0  0   0   0  0   0  -1  1
3×9 Array{Int8,2}:
 1  -1  1  0   0  0  0   0  0
 0   0  0  1  -1  1  0   0  0
 0   0  0  0   0  0  1  -1  1
```
"""
function uniteMultipleModels(models::Array{CAGD.Model,1})::CAGD.Model

    model = deepcopy(models[1])

    for i = 2 : length(models)
        CAGD.uniteModels!(model, models[i])
    end

    return model
end

"""
    mergeMultipleModels(models::Array{CAGD.Model,1}[; err=1e-6[, signed_merge=false]])::CAGD.Model

Unite `models` and performs CCE algorithm of the result with atol `err`.

If `signed_merge` is set true, `models` are expected signed and it returns a signed model.

See also:
 - [`CAGD.uniteMultipleModels`](@ref) for model unions;
 - [`CAGD.mergeModelVertices`](@ref) for CEE algorithm.
---
# Examples
```jldoctest
julia> models = [
    CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ]);

    CAGD.Model([0.0 0.0 -1.0; 0.0 -1.0 0.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ]);
    
    CAGD.Model([0.0 -1.0 -1.0; 1.0 0.0 1.0], [
        SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
        SparseArrays.sparse(Int8[1 -1 1])
    ])
];

julia> m = CAGD.mergeMultipleModels(models, signed_merge = true)
2D model with:
 - 6 points
 - 9 1-cells
 - 3 2-cells
2×6 Array{Float64,2}:
 0.0  0.0  1.0   0.0  -1.0  -1.0
 0.0  1.0  0.0  -1.0   0.0   1.0
9×6 Array{Int8,2}:
 -1   0  0   0   1  0
 -1   0  0   1   0  0
 -1   0  1   0   0  0
 -1   1  0   0   0  0
  0  -1  0   0   0  1
  0  -1  0   0   1  0
  0  -1  1   0   0  0
  0   0  0  -1   1  0
  0   0  0   0  -1  1
3×9 Array{Int8,2}:
 1  -1  0   0  0   0   0  -1   0
 0   0  1  -1  0   0  -1   0   0
 0   0  0   0  1  -1   0   0  -1

julia> Matrix(tgw(m, 2)[1])
4×9 Array{Int8,2}:
  1  -1   0   0   0   0   0  -1   0
  0   0   1  -1   0   0  -1   0   0
  0   0   0   0   1  -1   0   0  -1
 -1   0   0   1   0   1   0   0   0
  0   1  -1   0  -1   0   1   1   1
```
"""
function mergeMultipleModels(models::Array{CAGD.Model,1}; err=1e-6, signed_merge=false)::CAGD.Model
    return CAGD.mergeModelVertices(CAGD.uniteMultipleModels(models), err=err, signed_merge=signed_merge)
end

#-------------------------------------------------------------------------------
#   MODEL EXPORT
#-------------------------------------------------------------------------------

"""
    jlexportModel(model::CAGD.Model, filename::String[, modelname::String = "model"])::Nothing

Exports `model` by creating a file named `filename` that, if executed,
reloads the model with identifier `modelname` (wich is model, by default).
---
# Examples
```jldoctest
julia> modelOut = CAGD.Model([0.0 0.0 1.0; 0.0 1.0 0.0], [
    SparseArrays.sparse(Int8[-1 1 0; -1 0 1; 0 -1 1]),
    SparseArrays.sparse(Int8[1 -1 1])
]);

julia> CAGD.jlexportModel(modelOut, "test.jl", "modelIn")

julia> include("./test.jl");

julia> modelIn
2D model with:
 - 3 points
 - 3 1-cells
 - 1 2-cells
2×3 Array{Float64,2}:
 0.0  0.0  1.0
 0.0  1.0  0.0
3×3 Array{Int8,2}:
 -1   1  0
 -1   0  1
  0  -1  1
1×3 Array{Int8,2}:
 1  -1  1

julia> rm("./testModelExport.jl")
```
"""
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
