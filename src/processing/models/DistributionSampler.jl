struct DistributionSampler{N}
    shape_prior::ShapePrior{N};
    composition_smooth::Float64
end