@views maximize_noise_distributions(edge_lengths::Vector{Float64}, assignment_probs::Matrix{Float64}) =
    (Normal(wmean_std(edge_lengths, assignment_probs[:,i])...) for i in 1:2)

function expect_noise_probabilities!(
        assignment_probs::Matrix{Float64}, d1::Normal, d2::Normal, edge_lengths::Vector{Float64}, adj_list::AdjList;
        updating_ids::Vector{Int}, min_confidence::Union{Vector{Float64}, Nothing}=nothing
    )
    norm_denses = (pdf.(d1, edge_lengths), pdf.(d2, edge_lengths));
    n1 = sum(assignment_probs[:, 1])
    n2 = size(assignment_probs, 1) - n1

    for i in updating_ids
        cur_points = adj_list.ids[i]
        cur_weights = adj_list.weights[i]

        c_d1 = c_d2 = 0.0
        for j in 1:length(cur_points)
            ap = assignment_probs[cur_points[j], 1]
            c_d1 += cur_weights[j] * ap
            c_d2 += cur_weights[j] * (1. - ap)
        end

        d1 = n1 * exp(c_d1) * norm_denses[1][i]
        d2 = n2 * exp(c_d2) * norm_denses[2][i]
        p1 = d1 / fmax(d1 + d2, 1e-20)
        if min_confidence !== nothing
            p1 = min_confidence[i] + p1 * (1. - min_confidence[i])
        end

        assignment_probs[i, 1] = p1
        assignment_probs[i, 2] = 1. - p1
    end
end

function fit_noise_probabilities(
        edge_lengths::Vector{Float64}, adj_list::AdjList;
        min_confidence::Union{Vector{Float64}, Nothing}=nothing, max_iters::Int=10000, tol::Float64=0.005,
        verbose::Bool=false, progress::Union{Progress, Nothing}=nothing
    )
    # Initialization
    init_means = quantile(edge_lengths, (0.1, 0.9))
    init_std = (init_means[2] - init_means[1]) / 4.0
    d1, d2 = Normal(init_means[1], init_std), Normal(init_means[2], init_std)

    assignment_probs = hcat(pdf.(d1, edge_lengths), pdf.(d2, edge_lengths));

    ## Hard-assign outliers to the noise class to avoid numerical problems with super-small probabilities
    outlier_mask = (edge_lengths .> (init_means[2] + 3 * init_std))
    updating_ids = findall(.!outlier_mask)
    assignment_probs[outlier_mask, 2] .= 1.0
    assignment_probs ./= sum(assignment_probs, dims=2)
    max_diffs = Float64[]

    if verbose && progress === nothing
        progress = ProgressUnknown(0.3)
    end

    # EM iterations
    n_iters = max_iters
    for i in 1:max_iters
        expect_noise_probabilities!(
            assignment_probs, d1, d2, edge_lengths, adj_list; updating_ids, min_confidence
        )
        d1n, d2n = maximize_noise_distributions(edge_lengths, assignment_probs)

        ## Estimate parameter differences as convergence criteria
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
            n_iters = i
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

    if min_confidence !== nothing
        assignment_probs[:, 1] .= min_confidence .+ assignment_probs[:, 1] .* (1. .- min_confidence)
        assignment_probs[:, 2] .= 1.0 .- assignment_probs[:, 1]
    end

    assignment = Int.(assignment_probs[:, 2] .> 0.5) .+ 1

    return assignment_probs, assignment, (d1, d2), max_diffs
end

function estimate_confidence(df_spatial::DataFrame, prior_assignment::Union{Vector{Int}, Nothing}=nothing; nn_id::Int, prior_confidence::Float64=0.5)
    pos_data = position_data(df_spatial);
    mean_dists = getindex.(knn_parallel(KDTree(pos_data), pos_data, nn_id + 1; sorted=true)[2], nn_id + 1)

    min_confidence = (prior_assignment !== nothing) ? (prior_confidence^2 .* (prior_assignment .> 0)) : nothing

    adj_list = build_molecule_graph(df_spatial, filter=false); # TODO: can be optimized as we already have kNNs
    probs, (d1, d2) = fit_noise_probabilities(mean_dists, adj_list; min_confidence)[[1, 3]]

    return mean_dists, probs[:, 1], d1, d2
end

function append_confidence!(df_spatial::DataFrame, prior_assignment::Union{Vector{Int}, Symbol, Nothing}=:auto; nn_id::Int, prior_confidence::Float64=0.5)
    if prior_assignment isa Symbol
        if prior_assignment == :auto
            prior_assignment = (:prior_segmentation in propertynames(df_spatial)) ? df_spatial.prior_segmentation : nothing
        else
            prior_assignment = df_spatial[!, prior_assignment]
        end
    end

    mean_dists, df_spatial[!,:confidence], d1, d2 = estimate_confidence(df_spatial, prior_assignment, nn_id=nn_id, prior_confidence=prior_confidence)
    return mean_dists, df_spatial[!,:confidence], d1, d2
end