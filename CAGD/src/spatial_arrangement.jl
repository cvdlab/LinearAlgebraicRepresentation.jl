function spatial_arrangement(model; epsilon = 1e-6) begin
    split_model = CAGD.pairwise_decomposition(model)
    congr_model = CAGD.mergeModelVertices(split_model, err=epsilon, signed=true)
    giftw_model =
end
