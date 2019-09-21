import Random
import Plots

using Statistics

center_dists(pos_data::Matrix{Float64}, ids_per_clust::Array{Vector{Int}, 1}) =
    sum(diff(hcat([mean(pos_data[:, ids], dims=2) for ids in ids_per_clust]...), dims=2) .^ 2)

function permutation_mixing_score(cell_pos_data::Matrix{Float64}, cell_mol_clusts::Vector{Int}; n_permutations::Int=1000, threshold::Float64=0.05)
    ids_per_clust = split_ids(cell_mol_clusts)
    ids_per_clust = ids_per_clust[sortperm(length.(ids_per_clust), rev=true)]

    if length(ids_per_clust) == 1
        return 0.0
    end

    n_wrong_mols = 0
    for wrong_clust_ids in ids_per_clust[2:end]
        if length(wrong_clust_ids) < 2
            continue
        end

        obs_dist = center_dists(cell_pos_data, ids_per_clust)

        dists = Float64[]
        for i in 1:n_permutations
            s_ids = split_ids(Random.shuffle(cell_mol_clusts))
            push!(dists, center_dists(cell_pos_data, s_ids[sortperm(length.(s_ids), rev=true)[1:2]]))
        end

        if mean(dists .>= obs_dist) < threshold
            n_wrong_mols += length(wrong_clust_ids)
        end
    end

    return n_wrong_mols / length(cell_mol_clusts)
end

function mixing_score_per_cell(position_data::Matrix{Float64}, cell_assignment::Vector{Int}, clust_per_mol::Vector{Int}; min_mols_per_cell::Int)
    ids_per_cell = split_ids(cell_assignment .+ 1)[2:end]
    real_cell_ids = findall(length.(ids_per_cell) .>= min_mols_per_cell)
    wrong_frac_per_cell = [permutation_mixing_score(position_data[:, ids], clust_per_mol[ids]) for ids in ids_per_cell[real_cell_ids]];

    return wrong_frac_per_cell, length.(ids_per_cell[real_cell_ids])
end