using ProgressMeter
using DataFrames
using Statistics
using StatsBase
using Random

function maximize_mols!(cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, assignment::Vector{Int})
    cell_type_exprs .= 0.0

    for i in 1:length(assignment)
        cell_type_exprs[assignment[i], genes[i]] += 1.0
    end

    cell_type_exprs ./= sum(cell_type_exprs, dims=2)
end

function expect_mols!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.001,
        mrf_prior_weight::Float64=1.0) # Add confidence?
    denses = zeros(size(cell_type_exprs, 1))
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
                denses[assignment[cur_points[j]]] *= exp(cur_weights[j] * mrf_prior_weight)
            end

            assignment[i] = fsample(denses)
        end
    end
end

# TODO: rename
function optimize_mols(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        k::Int, n_iters::Int=100, history_depth::Int=50, new_prob::Float64=0.001, mrf_prior_weight::Float64=1.0)
    cell_type_exprs = copy(hcat([rand(maximum(genes)) for i in 1:k]...)');
    assignment = rand(1:k, length(genes));
    assignment_history = Vector{Int}[]
    maximize_mols!(cell_type_exprs, genes, assignment)
    @showprogress for i in 1:n_iters
        expect_mols!(assignment, cell_type_exprs, genes, adjacent_points, adjacent_weights, new_prob=new_prob,
            mrf_prior_weight=mrf_prior_weight)
        maximize_mols!(cell_type_exprs, genes, assignment)
        if (n_iters - i) >= history_depth
            push!(assignment_history, deepcopy(assignment))
        end
    end

    assignment_cons = vec(mapslices(mode, hcat(assignment_history...), dims=2));

    return cell_type_exprs, assignment_cons, assignment_history
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