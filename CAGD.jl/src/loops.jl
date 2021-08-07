"""
    getExternalCycle(model::CAGD.Model, dim::Int, cells::Array{Int,1})::Int

Evaluates wich `dim`-coboundary between `cells` is the external one.

The function expects `cells` as the indices of ALL the `dim`-coboundary
within the same component and evaluate wich of those is the external one.

It returns the index of the external cell within the `cells`.
"""
function getExternalCycle(model::CAGD.Model, dim::Int, cells::Array{Int,1})::Int
    
    return 1
end

"""
    buildContainmentGraph(
        model::CAGD.Model,
        dim::Int,
        externals::Array{Int,1},
        internals::Array{Array{Int,1},1}
    )::SparseMatrixCSC{Int, Int}

Evaluates the strict containment relation between `dim`-cells.

The method expects `externals[i]` and `internals[i]` as indices respectivelly
of the external and internal cycles of the `i`-th component out of `ne`.

Returns a `ne × ne` sparse `graph` where `graph[i, j] = k` if `i`-th component
includes `j`-th one within the `k`-th cycle (internal cycle of `i`-th component)
"""
function buildContainmentGraph(model, dim, externals, internals)
    ne = length(externals)

    # Bounding box for all the `dim`-cells
    bboxes = CAGD.getModelBoundingBoxes(model, dim)

    # graph[i, j] = 1 if bounding box of external[i] contains the bounding box of external[j]
    graph = SparseArrays.sparse([], [], Bool.([]), ne, ne)
    for c = 1 : ne  for r = 1 : ne
        if CAGD.isBoundingBoxIncluded(bboxes[externals[c]], bboxes[externals[r]])
            graph[r, c] = true
        end
    end  end

    # define a point per contained cell to be tested out
    pts = zeros(3, ne)
    containers, contained, _ = SparseArrays.findnz(graph)
    for c in unique(contained)
        pts[:, c] = CAGD.getModelCellVertices(model, dim, externals[c])[1]
    end
    
    # graph[r, c] = 1 if a point of c-th (external) is in r-th cell (external)
    for (r, c) in zip(containers, contained)
        if CAGD.pointInCell(model, dim, pts[c], externals[r], signed = false) >= 0
            SparseArrays.dropstored!(graph, r, c)
        end
    end
    
    # Transitive Reduction of graph
    containers, contained, _ = SparseArrays.findnz(graph)
    for (out, mid) in zip(containers, contained)
        for inn in SparseArrays.nzind(graph[mid, :])
            SparseArrays.dropstored!(graph, out, inn);
        end
    end

    # graph[r, c] = k if internal[r][k] contains external[c]
    containers, contained, _ = SparseArrays.findnz(graph)
    for (container, contained) in zip(containers, contained)
        for doContain in internals[container]
            if CAGD.isBoundingBoxIncluded(bboxes[externals[contained]], bboxes[doContain]) &
                    (CAGD.pointInCell(model, dim, pts[contained], doContain, signed = false) >= 0)
                graph[container, contained] = doContain
                break
            end
        end
    end
    
    return graph
end


"""
    assembleModel!(
        model::CAGD.Model,
        dim::Int,
        externals::Array{Int,1},
        containment_graph::SparseMatrixCSC{Int, Int},
        doExternal = false
    )::Nothing

Performs the adjoining of the outer cycles of inner components to the outer ones

Assemble on `model` the `ne` `dim`-cycles `externals` on inner cycle specified
by `contained_graph`. The latter is expected to be a `ne × ne` `SparseMatrixCSC`
where `containment_graph[i, j] = k` if `i`-th component contains the `j`-th
component within cycle `k`.

If `doExternal` is set to true, then the outer cycle is also adjoint.
"""
function assembleModel!(model, dim, externals, containment_graph, doExternal = false)
    containers, contained, containing = SparseArrays.findnz(containment_graph)

    # Adjoining outer cell of inner components to outer components
    for (out, inn) in zip(containing, contained)
        model.T[dim][out, :] += model.T[dim][externals[inn], :]
    end

    # Retrieve the external cycle (if needed)
    if doExternal
        external = SparseArrays.spzeros(Int8, size(model, dim, 2))
        for out in setdiff(1 : length(externals), unique(containers))
            external += model.T[dim][externals[out], :]
        end
        CAGD.addModelCells!(model, dim, external)
    end
    
    # Removing outer cells of components
    CAGD.deleteModelCells!(model, dim, externals)
    return
end