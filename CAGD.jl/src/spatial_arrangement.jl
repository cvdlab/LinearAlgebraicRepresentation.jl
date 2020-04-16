function spatial_arrangement(model; atol = 1e-6, outerCell = false)
    split_model = CAGD.pairwise_decomposition(model)

    congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)

    #bico_models = CAGD.split_biconnected_models(congr_model)
    gift_model, bicon_comps = CAGD.tgw(congr_model, 3, atol=atol)

    externals = [bc[CAGD.getExternalCoboundary(gift_model, 3, bc)] for bc in bicon_comps]
    internals = map(i -> setdiff(bicon_comps[i], externals[i]), i = 1 : length(externals))
    containment_graph = CAGD.buildContainmentGraph(gift_model, 3, externals, internals)
    
    # ?? REARRANGE OUTER CELLS

    assembled_model = CAGD.assembleModel(model, 3, externals, containment_graph, outerCell)

end
