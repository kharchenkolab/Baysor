@views maximize_noise_distributions(edge_lengths::Vector{Float64}, assignment_probs::Matrix{Float64}; updating_ids::T2 where T2 <: AbstractArray{Int} = 1:length(edge_lengths)) =
    (Normal(wmean_std(edge_lengths, assignment_probs[:,i], non_zero_ids=updating_ids)...) for i in 1:2)

function expect_noise_probabilities!(assignment_probs::Matrix{Float64}, d1::Normal, d2::Normal, edge_lengths::Vector{Float64},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; updating_ids::T where T <: AbstractVector{Int} = 1:length(edge_lengths))
    norm_denses = (pdf.(d1, edge_lengths), pdf.(d2, edge_lengths));
    n1 = sum(assignment_probs[:, 1])
    n2 = size(assignment_probs, 1) - n1

    dists = (d1, d2)
    for i in updating_ids
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        c_d = 0.0
        sum_weights = 0.0
        for j in 1:length(cur_points)
            c_d += cur_weights[j] * assignment_probs[cur_points[j], 1]
            sum_weights += cur_weights[j]
        end

        d1 = n1 * exp(c_d) * norm_denses[1][i]
        d2 = n2 * exp(sum_weights - c_d) * norm_denses[2][i]
        ds = fmax(d1 + d2, 1e-20)

        assignment_probs[i, 1] = d1 / ds
        assignment_probs[i, 2] = d2 / ds
    end
end

function fit_noise_probabilities(edge_lengths::Vector{Float64}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        max_iters::Int=1000, tol::Float64=0.005, verbose::Bool=false, progress::Union{Progress, Nothing}=nothing)
    init_means = quantile(edge_lengths, (0.1, 0.9))
    init_std = (init_means[2] - init_means[1]) / 4.0
    d1, d2 = Normal(init_means[1], init_std), Normal(init_means[2], init_std)

    assignment_probs = hcat(pdf.(d1, edge_lengths), pdf.(d2, edge_lengths));
    outlier_mask = (edge_lengths .> (init_means[2] + 3 * init_std))
    updating_ids = findall(.!outlier_mask)
    assignment_probs[outlier_mask, 2] .= 1.0
    assignment_probs ./= sum(assignment_probs, dims=2)
    max_diffs = Float64[]

    if verbose && progress === nothing
        progress = ProgressUnknown(0.3)
    end

    edge_ord = sortperm(edge_lengths)
    n_iters = 0
    for i in 1:max_iters
        n_iters = i
        expect_noise_probabilities!(assignment_probs, d1, d2, edge_lengths, adjacent_points, adjacent_weights, updating_ids=updating_ids)
        d1n, d2n = maximize_noise_distributions(edge_lengths, assignment_probs, updating_ids=updating_ids)

        param_diff = max(abs(d1n.μ - d1.μ) / d1.μ, abs(d2n.μ - d2.μ) / d2.μ, abs(d1n.σ - d1.σ) / d1.σ, abs(d2n.σ - d2.σ) / d2.σ)
        push!(max_diffs, param_diff)
        d1, d2 = d1n, d2n

        if verbose
            next!(progress)
        end

        if max_diffs[end] < tol
            if verbose
                finish!(progress)
            end
            break
        end
    end

    if verbose
        println("Algorithm stopped after $n_iters iterations. Error: $(round(max_diffs[end], sigdigits=3)). Converged: $(max_diffs[end] <= tol).")
    end

    if d1.μ > d2.μ
        d2, d1 = d1, d2
    end

    n1 = sum(assignment_probs[:, 1])
    n2 = size(assignment_probs, 1) - n1
    assignment_probs = hcat(n1 .* pdf.(d1, edge_lengths), n2 .* pdf.(d2, edge_lengths))
    assignment_probs[outlier_mask, 2] .= 1.0
    assignment_probs ./= sum(assignment_probs, dims=2)

    assignment = Int.(assignment_probs[:, 2] .> 0.5) .+ 1

    return assignment_probs, assignment, (d1, d2), max_diffs
end

"""
Parameters:
- min_noise_level: minimal fraction of molecules, which are hard assigned to the noise class
- max_noise_level: maximal fraction of molecules, which are hard assigned to the noise class
- min_real_level: minimal fraction of molecules, which are hard assigned to the real class
"""
function append_confidence!(df_spatial::DataFrame, segmentation_mask::Union{BitArray{2}, Nothing}=nothing; nn_id::Int,
        border_quantiles::Tuple{Float64, Float64}=(0.3, 0.975), seg_quantile_left::Float64=0.9, seg_quantile_mid::Float64=0.99)
    pos_data = position_data(df_spatial);
    mean_dists = getindex.(knn(KDTree(pos_data), pos_data, nn_id + 1, true)[2], nn_id + 1)

    if segmentation_mask !== nothing
        border_left, border_right = quantile(mean_dists, border_quantiles)
        mol_in_nuclei_mask = segmentation_mask[CartesianIndex.(round.(Int, pos_data[2,:]), round.(Int, pos_data[1,:]))]
        sb_left, sb_mid = quantile(mean_dists[mol_in_nuclei_mask], [seg_quantile_left, seg_quantile_mid])
        border_left = max(sb_left, border_left) # TODO: maybe estimate with DAPI will be more robust, so should remove min/max
        border_right = min(sb_left + 2 * (sb_mid - sb_left), border_right)
        df_spatial[!,:confidence] = interpolate_linear.(mean_dists, border_left, border_right);
        return df_spatial.confidence
    end

    adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false); # TODO: can be optimized as we already have kNNs
    probs, (d1, d2) = fit_noise_probabilities(mean_dists, adjacent_points, adjacent_weights)[[1, 3]]
    df_spatial[!,:confidence] = probs[:, 1];
end