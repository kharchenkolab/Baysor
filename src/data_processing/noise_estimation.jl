function maximize_noise_distributions(edge_lengths::Vector{Float64}, assignment_probs::Matrix{Float64}; ord::Vector{Int}=sortperm(edge_lengths))
    m1, m2 = [wmedian(edge_lengths, assignment_probs[i,:], ord=ord) for i in 1:2]
    s1, s2 = [wmad(edge_lengths, assignment_probs[i,:], (m1, m2)[i]) for i in 1:2]
    return Normal(m1, s1), Normal(m2, s2)
end

function expect_noise_probabilities!(assignment_probs::Matrix{Float64}, d1::Normal, d2::Normal, edge_lengths::Vector{Float64},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}})
    norm_denses = copy(hcat(pdf.(d1, edge_lengths), pdf.(d2, edge_lengths))');
    dists = (d1, d2)
    for i in 1:length(edge_lengths)
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        assignment_probs[:, i] .= 0.0
        dense_sum = 0.0
        for ri in 1:2
            c_d = 0.0
            for j in 1:length(cur_points)
                c_d += cur_weights[j] * assignment_probs[ri, cur_points[j]]
            end

            assignment_probs[ri, i] = exp(c_d) * pdf(dists[ri], edge_lengths[i])
            dense_sum += assignment_probs[ri, i]
        end

        if dense_sum > 1e-20
            assignment_probs[1, i] /= dense_sum
            assignment_probs[2, i] /= dense_sum
        end
    end
end

function fit_noise_probabilities(edge_lengths::Vector{Float64}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        max_iters::Int=1000, tol::Float64=0.01, verbose::Bool=false, progress::Union{Progress, Nothing}=nothing)
    init_means = quantile(edge_lengths, (0.1, 0.9))
    init_std = (init_means[2] - init_means[1]) / 4.0
    d1, d2 = Normal(init_means[1], init_std), Normal(init_means[2], init_std)

    assignment_probs = copy(hcat(pdf.(d1, edge_lengths), pdf.(d2, edge_lengths))');
    assignment_probs ./= sum(assignment_probs, dims=1)
    assignment_probs_prev = deepcopy(assignment_probs)
    max_diffs = Float64[]

    if verbose && progress === nothing
        progress = ProgressUnknown(0.3)
    end

    edge_ord = sortperm(edge_lengths)
    n_iters = 0
    for i in 1:max_iters
        n_iters = i
        assignment_probs_prev .= assignment_probs
        expect_noise_probabilities!(assignment_probs, d1, d2, edge_lengths, adjacent_points, adjacent_weights)
        d1, d2 = maximize_noise_distributions(edge_lengths, assignment_probs, ord=edge_ord)

        push!(max_diffs, Baysor.estimate_difference_l0(assignment_probs, assignment_probs_prev))
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

    assignment = Int.(assignment_probs[2, :] .> assignment_probs[1, :]) .+ 1

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
    confidence = d1.μ > d2.μ ? probs[2, :] : probs[1, :]
    df_spatial[!,:confidence] = confidence;
end