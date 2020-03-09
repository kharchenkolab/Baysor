using ProgressMeter
using DataFrames
using Statistics
using StatsBase
using Random

function maximize_mols!(cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, assignment::Vector{Int})
    cell_type_exprs[unique(assignment), :] .= 0.0

    for i in 1:length(assignment)
        cell_type_exprs[assignment[i], genes[i]] += 1.0
    end
end

function expect_mols!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.001) # Add confidence?
    denses = zeros(size(cell_type_exprs, 1))
    cell_type_exprs = (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1) # It's not sum to 1, but it simulates sparse pseudocounts
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        denses .= 0.0

        if rand() < new_prob # allowing to link new component improves global convergence
            assignment[i] = fsample(cell_type_exprs[:, gene])
        else
            for j in 1:length(cur_points)
                ca = assignment[cur_points[j]]
                denses[ca] = cell_type_exprs[ca, gene]
            end

            for j in 1:length(cur_weights)
                denses[assignment[cur_points[j]]] *= cur_weights[j]
            end

            assignment[i] = fsample(denses)
        end
    end
end

function filter_correlated_clusters!(cell_type_exprs::Matrix{Float64}, assignment::Vector{Int}; correlation_threshold::Float64=0.95)
    cors = cor(cell_type_exprs');
    cors[diagind(cors)] .= 0
    triu!(cors);
    max_cor_ids = vec(mapslices(findmax, cors, dims=1));

    was_filtering = false
    for i1 in 1:length(max_cor_ids)
        c, i2 = max_cor_ids[i1]
        if c < correlation_threshold
            continue
        end

        was_filtering = true
        cell_type_exprs[i1, :] .= 0
        assignment[assignment .== i1] .= i2
    end

    return was_filtering
end

function remove_unused_clusters!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int}; min_mols_per_type)
    n_mols_per_type = count_array(assignment)
    real_type_ids = findall(n_mols_per_type .>= min_mols_per_type)
    if length(real_type_ids) == length(n_mols_per_type)
        return cell_type_exprs
    end

    id_map = zeros(Int, size(cell_type_exprs, 1))
    id_map[real_type_ids] .= 1:length(real_type_ids)

    cell_type_exprs = cell_type_exprs[real_type_ids,:]

    cell_type_exprs_norm = cell_type_exprs ./ sum(cell_type_exprs, dims=2)
    for i in 1:length(assignment)
        if n_mols_per_type[assignment[i]] < min_mols_per_type
            assignment[i] = findmax(cell_type_exprs_norm[:, genes[i]])[2]
        else
            assignment[i] = id_map[assignment[i]]
        end
    end

    return cell_type_exprs
end

# TODO: rename
function optimize_mols(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        k::Int, n_iters::Int=1000, history_depth::Int=200, new_prob::Float64=0.05, mrf_prior_weight::Float64=1.0,
        min_mols_per_type::Int = round(Int, 0.05 * length(genes) / k), correlation_threshold::Float64=0.95,
        verbose::Bool=true, progress::Union{Progress, Nothing}=nothing)
    cell_type_exprs = zeros(k, maximum(genes));
    assignment = rand(1:k, length(genes));
    assignment_history = Vector{Int}[]
    maximize_mols!(cell_type_exprs, genes, assignment)

    adjacent_weights = [exp.(aw .* mrf_prior_weight) for aw in adjacent_weights]

    if verbose && progress === nothing
        progress = Progress(n_iters, 0.3)
    end

    for i in 1:n_iters
        expect_mols!(assignment, cell_type_exprs, genes, adjacent_points, adjacent_weights, new_prob=new_prob)
        maximize_mols!(cell_type_exprs, genes, assignment)
        if (n_iters - i) == (history_depth - 1)
            cell_type_exprs = remove_unused_clusters!(assignment, cell_type_exprs, genes; min_mols_per_type=min_mols_per_type)
            new_prob = 0.0
        end

        if (n_iters - i) < history_depth
            push!(assignment_history, deepcopy(assignment))
        else
            if filter_correlated_clusters!(cell_type_exprs, assignment, correlation_threshold=correlation_threshold)
                maximize_mols!(cell_type_exprs, genes, assignment)
            end
        end

        if verbose
            next!(progress)
        end
    end

    assignment_cons = vec(mapslices(mode, hcat(assignment_history...), dims=2));

    return cell_type_exprs ./ sum(cell_type_exprs, dims=2), assignment_cons, assignment_history
end

function build_molecule_graph(df_spatial::DataFrame; kwargs...)
    edge_list, adjacent_dists = adjacency_list(df_spatial; kwargs...);

    real_edge_length = quantile(adjacent_dists, 0.3);
    adjacent_weights = real_edge_length ./ max.(adjacent_dists, real_edge_length);

    return convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1));
end

# function build_molecule_knn_graph(df_spatial::DataFrame)
#     edge_list, adjacent_dists = knn(KDTree(position_data(df_spatial)), position_data(df_spatial), 6, true)
#     edge_list = [el[2:end] for el in edge_list]
#     adjacent_dists = [ad[2:end] for ad in adjacent_dists]

#     real_edge_length = quantile(vcat(adjacent_dists...), 0.3);
#     adjacent_weights = [real_edge_length ./ max.(ad, real_edge_length) for ad in adjacent_dists];

#     return edge_list, adjacent_weights
# end