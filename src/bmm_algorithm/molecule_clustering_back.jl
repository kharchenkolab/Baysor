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

function generage_new_component(genes::Vector{Int}, adj_points::Vector{Int}; n_genes::Int)
    gene_probs = zeros(1, n_genes)
    @inbounds for p in adj_points
        gene_probs[1, genes[p]] += 1. / length(adj_points)
    end

    return gene_probs
end

function expect_mols!(assignment::Vector{Int}, n_mols_per_comp::Vector{Int}, gene_comp_probs::Vector{Float64}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.001, α::Float64=1.0) # Add confidence?
    denses = zeros(size(cell_type_exprs, 1))
    med_weight = median(vcat(adjacent_weights...))
    # quant_weight = quantile(vcat(adjacent_weights...), α)
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        sum_dens = sum(cur_weights)
        dens_thres = sum_dens * 1e-7
        denses .= 0.0
        new_c = 0

        if rand() < new_prob
            new_c = fsample(cell_type_exprs[:, gene])
        else
#             n_mols_per_comp[assignment[i]] -= 1
            @inbounds @simd for j in 1:length(cur_weights)
                denses[assignment[cur_points[j]]] += cur_weights[j]
            end

            @inbounds @simd for j in 1:length(denses)
                if denses[j] > dens_thres
                    # denses[j] *= log2(max(n_mols_per_comp[j], α)) * max(cell_type_exprs[j, gene], 1e-2) # TODO: better smoothing
                    # denses[j] *= max(n_mols_per_comp[j], α) * max(cell_type_exprs[j, gene], 1e-2) # TODO: better smoothing
                    denses[j] *= max(cell_type_exprs[j, gene], 1e-3) # TODO: better smoothing
                end
            end
            push!(denses, α * med_weight * gene_comp_probs[i] / length(genes))
            # push!(denses, quant_weight * gene_comp_probs[i])

            new_c = fsample(denses)
            if new_c == length(denses)
                cell_type_exprs = vcat(cell_type_exprs, generage_new_component(genes, cur_points, n_genes=size(cell_type_exprs, 2)))
                push!(n_mols_per_comp, 0)
            else
                pop!(denses)
            end
        end

        n_mols_per_comp[assignment[i]] -= 1
        assignment[i] = new_c
        n_mols_per_comp[new_c] += 1
    end

    return cell_type_exprs
end

# TODO: rename
function optimize_mols(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        k::Int, n_iters::Int=100, history_depth::Int=50, new_prob::Float64=0.001, α::Float64=1.0)
    cell_type_exprs = copy(hcat([rand(maximum(genes)) for i in 1:k]...)');
    assignment = rand(1:k, length(genes));
    assignment_history = Vector{Int}[]
    maximize_mols!(cell_type_exprs, genes, assignment)

    n_mols_per_comp = count_array(assignment)
    gene_comp_probs = [(sum(genes[adjacent_points[i]] .== genes[i]) + 1) / (length(adjacent_points[i]) + 1) for i in 1:length(genes)]

    @showprogress for i in 1:n_iters
        cell_type_exprs = expect_mols!(assignment, n_mols_per_comp, gene_comp_probs, cell_type_exprs, genes, adjacent_points, adjacent_weights,
            new_prob=new_prob, α=α)
        maximize_mols!(cell_type_exprs, genes, assignment)
        if (n_iters - i) >= history_depth
            push!(assignment_history, deepcopy(assignment))
        end

        count_array!(n_mols_per_comp, assignment)
        non_empty_ids = findall(n_mols_per_comp .> 0)
        if length(non_empty_ids) < length(n_mols_per_comp)
            id_map = zeros(Int, size(cell_type_exprs, 1))
            id_map[non_empty_ids] .= 1:length(non_empty_ids)

            n_mols_per_comp = n_mols_per_comp[non_empty_ids]
            cell_type_exprs = cell_type_exprs[non_empty_ids, :]

            assignment = id_map[assignment]
        end
    end

    # TODO: adjust in history
    # assignment_cons = vec(mapslices(mode, hcat(assignment_history...), dims=2));

    return cell_type_exprs, assignment, n_mols_per_comp, assignment_history
    # return cell_type_exprs, assignment_cons, n_mols_per_comp, assignment_history
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
