function spatial_arrangement(model; atol = 1e-6, outerCell = false)
    split_model = CAGD.pairwise_decomposition(model, atol = atol)

    congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)

    #bico_models = CAGD.split_biconnected_models(congr_model)
    FCcycles, bicon_comps = CAGD.tgw(congr_model, 3, atol=atol)

    externals = [bc[CAGD.getExternalBoundary(congr_model, 3, FCcycles, bc)] for bc in bicon_comps]
    internals = [setdiff(bicon_comps[i], externals[i]) for i = 1 : length(externals)]
    # containment_graph = CAGD.buildContainmentGraph(gift_model, 3, externals, internals)
    
    # ?? REARRANGE OUTER CELLS

    # assembled_model = CAGD.assembleModel(model, 3, externals, containment_graph, outerCell)

    # WILL BE REMOVED
    # Sums all of the component external in a single component
    assembled_model = deepcopy(congr_model)

    FC = FCcycles[:, internals...]
    if outerCell
        JK = ones(Int8, length(externals))
        FC = [FCcycles * SparseArrays.sparse(externals, JK, JK, size(FCcycles, 2), 1) FC]
    end
    CAGD.addModelCells!(assembled_model, 3, convert(Lar.ChainOp, FC'))
    return assembled_model

end
