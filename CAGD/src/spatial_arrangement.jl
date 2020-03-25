function spatial_arrangement(model; epsilon = 1e-6)
    split_model = CAGD.pairwise_decomposition(model)
    congr_model = CAGD.mergeModelVertices(split_model, err=epsilon, signed=true)

    #bico_models = CAGD.split_biconnected_models(congr_model)

    gift_models = CAGD.TGW(congr_model, 3)
end
