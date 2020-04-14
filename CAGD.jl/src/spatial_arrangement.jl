function spatial_arrangement(model; atol = 1e-6)
    split_model = CAGD.pairwise_decomposition(model)
    congr_model = CAGD.mergeModelVertices(split_model, err=atol, signed_merge=true)
    
    #bico_models = CAGD.split_biconnected_models(congr_model)

    gift_model, bicon_comps = CAGD.tgw(congr_model, 3, atol=atol)
end
