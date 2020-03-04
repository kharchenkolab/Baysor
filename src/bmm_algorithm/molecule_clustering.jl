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

# function expect_mols!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
#         adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}) # Add confidence?
#     component_weights = Dict{Int, Float64}();
#     adj_classes = Int[]
#     adj_weights = Float64[]
#     denses = Float64[]
#     for i in 1:length(genes)
#         gene = genes[i]
#         adjacent_component_weights!(adj_weights, adj_classes, component_weights, assignment, adjacent_points[i], adjacent_weights[i])

#         empty!(denses)
#         for j in eachindex(adj_weights)
#             push!(denses, adj_weights[j] * cell_type_exprs[adj_classes[j], gene])
#         end

#         assignment[i] = fsample(adj_classes, denses)
#     end
# end

function expect_mols!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.05) # Add confidence?
    denses = zeros(size(cell_type_exprs, 1))
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        sum_dens = sum(cur_weights)
        denses .= 0.0

        if rand() < new_prob
            assignment[i] = fsample(cell_type_exprs[:, gene])
        else
            @inbounds @simd for j in 1:length(cur_weights)
                denses[assignment[cur_points[j]]] += cur_weights[j]
            end
            denses .*= cell_type_exprs[:, gene]

            assignment[i] = fsample(denses)
        end
    end
end

# TODO: rename
function optimize_mols(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        k::Int, n_iters::Int=100, history_depth::Int=50, new_prob::Float64=0.001)
    cell_type_exprs = copy(hcat([rand(maximum(genes)) for i in 1:k]...)');
    assignment = rand(1:k, length(genes));
    assignment_history = Vector{Int}[]
    maximize_mols!(cell_type_exprs, genes, assignment)
    @showprogress for i in 1:n_iters
        expect_mols!(assignment, cell_type_exprs, genes, adjacent_points, adjacent_weights, new_prob=new_prob)
        maximize_mols!(cell_type_exprs, genes, assignment)
        if (n_iters - i) >= history_depth
            push!(assignment_history, deepcopy(assignment))
        end
    end

    assignment_cons = vec(mapslices(mode, hcat(assignment_history...), dims=2));

    return cell_type_exprs, assignment_cons, assignment_history
end

function build_molecule_graph(df_spatial::DataFrame)
    edge_list, adjacent_dists = adjacency_list(df_spatial);
    real_edge_length = quantile(adjacent_dists, 0.3);
    adjacent_weights = 1 ./ max.(adjacent_dists, real_edge_length);

    return convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1));
end

# function reorder_clusters(cell_type_exprs::Matrix{Float64}, assignment::Vector{Int})
#     cl_order = sortperm(count_array(assignment), rev=true);
#     return cell_type_exprs[cl_order,:], get.(Ref(Dict(Pair.(cl_order, 1:10))), assignment, 0)
# end

# reorder_clusters(assignment::Vector{Int}) =
#     get.(Ref(Dict(Pair.(sortperm(count_array(assignment), rev=true), 1:10))), assignment, 0)
