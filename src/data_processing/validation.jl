import Random
import Plots

using Statistics

center_dists(pos_data::Matrix{Float64}, ids_per_clust1::Vector{Int}, ids_per_clust2::Vector{Int}) =
    sum(diff(hcat([mean(pos_data[:, ids], dims=2) for ids in [ids_per_clust1, ids_per_clust2]]...), dims=2) .^ 2)

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

        obs_dist = center_dists(cell_pos_data, ids_per_clust[1], wrong_clust_ids)

        dists = Float64[]
        for i in 1:n_permutations
            labels_shuffled = Random.shuffle(cell_mol_clusts)
            push!(dists, center_dists(cell_pos_data,
                findall(labels_shuffled .== cell_mol_clusts[ids_per_clust[1][1]]),
                findall(labels_shuffled .== cell_mol_clusts[wrong_clust_ids[1]])))
        end

        if mean(dists .>= obs_dist) < threshold
            n_wrong_mols += length(wrong_clust_ids)
        end
    end

    return n_wrong_mols / length(cell_mol_clusts)
end

function mixing_score_per_cell(position_data::Matrix{Float64}, cell_assignment::Vector{Int}, clust_per_mol::Vector{Int}; min_mols_per_cell::Int, kwargs...)
    ids_per_cell = split_ids(cell_assignment .+ 1)[2:end]
    real_cell_ids = findall(length.(ids_per_cell) .>= min_mols_per_cell)
    wrong_frac_per_cell = [permutation_mixing_score(position_data[:, ids], clust_per_mol[ids]; kwargs...) for ids in ids_per_cell[real_cell_ids]];

    return wrong_frac_per_cell, length.(ids_per_cell[real_cell_ids])
end

function mixing_score_per_molecule(position_data::Matrix{Float64}, cell_assignment::Vector{Int}, clust_per_mol::Vector{Int}; min_mols_per_cell::Int, kwargs...)
    ids_per_cell = split_ids(cell_assignment .+ 1)[2:end]
    real_cell_ids = findall(length.(ids_per_cell) .>= min_mols_per_cell)
    wrong_frac_per_cell = [permutation_mixing_score(position_data[:, ids], clust_per_mol[ids]; kwargs...) for ids in ids_per_cell[real_cell_ids]];
    wrong_frac_per_mol = zeros(length(cell_assignment))
    for (frac, ids) in zip(wrong_frac_per_cell, ids_per_cell[real_cell_ids])
        wrong_frac_per_mol[ids] .= frac
    end

    return wrong_frac_per_mol
end

function plot_mixing_scores(cell_score_baysor::Vector{Float64}, n_mols_per_cell_baysor::Vector{Int},
                            cell_score_paper::Vector{Float64}, n_mols_per_cell_paper::Vector{Int},
                            cell_score_random::Union{Vector{Float64}, Nothing}, n_mols_per_cell_random::Union{Vector{Int}, Nothing}=nothing;
                            min_mols_per_cell::Int=0, layout=(2, 2), legend=:best, size1=(400, 300), size2=(800, 600), kwargs...)

    cell_score_baysor = cell_score_baysor[n_mols_per_cell_baysor .>= min_mols_per_cell]
    n_mols_per_cell_baysor = n_mols_per_cell_baysor[n_mols_per_cell_baysor .>= min_mols_per_cell]
    cell_score_paper = cell_score_paper[n_mols_per_cell_paper .>= min_mols_per_cell]
    n_mols_per_cell_paper = n_mols_per_cell_paper[n_mols_per_cell_paper .>= min_mols_per_cell]

    if cell_score_random !== nothing
        cell_score_random = cell_score_random[n_mols_per_cell_random .>= min_mols_per_cell]
        n_mols_per_cell_random = n_mols_per_cell_random[n_mols_per_cell_random .>= min_mols_per_cell]
    end

    p0 = Plots.scatter(n_mols_per_cell_baysor, cell_score_baysor, xlabel="Cell size", ylabel="Mixing score", label="Baysor",
        markersize=3, size=size1, legend=legend)
    p0 = Plots.scatter!(n_mols_per_cell_paper, cell_score_paper, label="Paper", markersize=3)

    if cell_score_random !== nothing
        p0 = Plots.scatter!(n_mols_per_cell_random, cell_score_random, label="Random", markersize=3)
    end

    t_xvals = sort(unique(n_mols_per_cell_paper))
    p1 = Plots.scatter(t_xvals, [mean((cell_score_paper .* n_mols_per_cell_paper)[n_mols_per_cell_paper .<= v]) for v in t_xvals], label="Paper",
        legend=legend, bg_legend=Colors.RGBA(1, 1, 1, 0.5), fg_legend=Colors.RGBA(0, 0, 0, 0.0), xlabel="Max cell size", ylabel="Mean number of 'wrong'\nmolecules per cell")

    t_xvals = sort(unique(n_mols_per_cell_baysor))
    p1 = Plots.scatter!(t_xvals, [mean((cell_score_baysor .* n_mols_per_cell_baysor)[n_mols_per_cell_baysor .<= v]) for v in t_xvals], label="Baysor")

    if cell_score_random !== nothing
        t_xvals = sort(unique(n_mols_per_cell_random))
        p1 = Plots.scatter!(t_xvals, [mean((cell_score_random .* n_mols_per_cell_random)[n_mols_per_cell_random .<= v]) for v in t_xvals], label="Random")
    end

    t_xvals = sort(unique(n_mols_per_cell_paper))
    p2 = Plots.scatter(t_xvals, [mean(cell_score_paper[n_mols_per_cell_paper .<= v]) for v in t_xvals], label="Paper",
        legend=:none, xlabel="Max cell size", ylabel="Contamination", ylim=(0.0, 1.0))

    t_xvals = sort(unique(n_mols_per_cell_baysor))
    p2 = Plots.scatter!(t_xvals, [mean(cell_score_baysor[n_mols_per_cell_baysor .<= v]) for v in t_xvals], label="Baysor")

    t_xvals = sort(unique(n_mols_per_cell_paper))
    p3 = Plots.scatter(t_xvals, [mean((cell_score_paper .* n_mols_per_cell_paper)[n_mols_per_cell_paper .>= v]) for v in t_xvals], label="Paper",
        legend=:none, xlabel="Min cell size", ylabel="Mean number of 'wrong'\nmolecules per cell")

    t_xvals = sort(unique(n_mols_per_cell_baysor))
    p3 = Plots.scatter!(t_xvals, [mean((cell_score_baysor .* n_mols_per_cell_baysor)[n_mols_per_cell_baysor .>= v]) for v in t_xvals], label="Baysor")

    if cell_score_random !== nothing
        t_xvals = sort(unique(n_mols_per_cell_random))
        p3 = Plots.scatter!(t_xvals, [mean((cell_score_random .* n_mols_per_cell_random)[n_mols_per_cell_random .>= v]) for v in t_xvals], label="Random")
    end

    t_xvals = sort(unique(n_mols_per_cell_paper))
    p4 = Plots.scatter(t_xvals, [mean(cell_score_paper[n_mols_per_cell_paper .>= v]) for v in t_xvals], label="Paper",
        legend=:none, xlabel="Min cell size", ylabel="Contamination", ylim=(0.0, 1.0))

    t_xvals = sort(unique(n_mols_per_cell_baysor))
    p4 = Plots.scatter!(t_xvals, [mean(cell_score_baysor[n_mols_per_cell_baysor .>= v]) for v in t_xvals], label="Baysor")

    if cell_score_random !== nothing
        t_xvals = sort(unique(n_mols_per_cell_random))
        p4 = Plots.scatter!(t_xvals, [mean(cell_score_random[n_mols_per_cell_random .>= v]) for v in t_xvals], label="Random")
    end

    return p0, Plots.plot(p1, p2, p3, p4; layout=layout, size=size2, kwargs...)
end

function assignment_summary_df(assignments::Pair{Symbol, Vector{Int}}...; min_molecules_per_cell::Int)
    assignments = Dict(assignments)
    return DataFrame(Dict(
        "name" => collect(keys(assignments)),
        "noise_frac" => [round(mean(x .== 0), sigdigits=3) for x in values(assignments)],
        "n_cells" => [sum(count_array(x .+ 1)[2:end] .>= min_molecules_per_cell) for x in values(assignments)]
    ))[:, [:name, :n_cells, :noise_frac]]
end

function convert_segmentation_to_counts(genes::Vector{Int}, cell_assignment::Vector{Int}; drop_empty_labels::Bool=false)
    if drop_empty_labels
        if minimum(cell_assignment) == 0
            cell_assignment = denserank(cell_assignment) .- 1
        else
            cell_assignment = denserank(cell_assignment)
        end
    end

    cm = zeros(Int, maximum(genes), maximum(cell_assignment))
    for i in 1:length(genes)
        if cell_assignment[i] == 0
            continue
        end
        cm[genes[i], cell_assignment[i]] += 1
    end

    return cm
end