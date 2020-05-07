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

function convert_segmentation_to_counts(genes::Vector{Int}, cell_assignment::Vector{Int}; drop_empty_labels::Bool=false, gene_names::Union{Vector{String}, Nothing}=nothing)
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

    if gene_names !== nothing
        cm = DataFrame(cm, [Symbol("$c") for c in 1:size(cm, 2)])
        cm[!, :gene] = gene_names
        cm = cm[:, vcat(end, 1:end-1)]
    end

    return cm
end

function plot_subset(df_spatial::DataFrame, dapi_arr::Matrix{<:Real}, (xs, xe), (ys, ye); polygons::Union{Bool, Vector{Matrix{Float64}}}=true, ms=2.0, alpha=0.2,
        grid_step::Float64=5.0, bandwidth::Float64=grid_step, cell_col::Symbol=:cell, dapi_alpha=0.9, polygon_line_width=2, noise::Bool=true, size_mult=1/3,
        plot_raw_dapi::Bool=true, color_col::Symbol=:color, annotation_col::Union{Symbol, Nothing}=nothing, build_panel::Bool=true, grid_alpha::Float64=0.5, ticks=false, kwargs...)
    df_subs = @where(df_spatial, :x .>= xs, :x .<= xe, :y .>= ys, :y .<= ye);

    if (typeof(polygons) == Bool)
        if polygons
            polygons = boundary_polygons(df_subs, df_subs[!, cell_col], grid_step=grid_step, min_molecules_per_cell=10, bandwidth=bandwidth)
        else
            polygons = Matrix{Float64}[]
        end
    end

    # xs, xe, ys, ye = round.(Int, [minimum(df_subs.x), maximum(df_subs.x), minimum(df_subs.y), maximum(df_subs.y)]);

    xticks_vals = range(0, xe-xs, length=5)
    yticks_vals = range(0, ye-ys, length=5)

    yticks = xticks = false
    if ticks
        xticks = (xticks_vals, ["$f" for f in round.(Int, range(xs, xe, length=5))])
        yticks = (yticks_vals, ["$f" for f in round.(Int, range(ys, ye, length=5))])
    end

    dapi_subs = dapi_arr[ys:ye, xs:xe]
    plot_size = ((xe-xs), ye-ys) .* size_mult
    plt1 = Plots.heatmap(maximum(dapi_subs) .- dapi_subs, color=:grayscale, colorbar=:none, size=plot_size,
                alpha=dapi_alpha, format=:png, legend=:none, xticks=xticks, yticks=yticks)

    is_noise = noise ? (df_subs[!, cell_col] .== 0) : nothing

    annotation = (annotation_col === nothing) ? nothing : df_subs[!, annotation_col]
    plot_cell_borders_polygons!(df_subs, polygons; color=df_subs[!, color_col], ms=ms, alpha=alpha, offset=(-xs, -ys),
        polygon_line_width=polygon_line_width, polygon_alpha=0.75, is_noise=is_noise, noise_kwargs=Dict(:ms => 1.0), annotation=annotation, kwargs...)

    Plots.vline!(xticks_vals, color="black", alpha=grid_alpha)
    Plots.hline!(yticks_vals, color="black", alpha=grid_alpha)

    if !plot_raw_dapi
        return plt1
    end

    # Possible colorschemes: tarn, diff, lime_grad, thermal
    plt2 = Plots.heatmap(dapi_subs, color=:diff, colorbar=:none, alpha=0.9, format=:png, xticks=xticks, yticks=yticks, legend=:none, size=plot_size)

    Plots.plot!([Plots.Shape(pg[:,1] .- xs, pg[:,2] .- ys) for pg in polygons],
        fill=(0, 0.0), linewidth=2.0, linecolor="black", alpha=0.4, label="", xlims=(0, (xe-xs)), ylims=(0, (ye-ys)));

    Plots.vline!(xticks_vals, color="black", alpha=grid_alpha)
    Plots.hline!(yticks_vals, color="black", alpha=grid_alpha)

    if !build_panel
        return plt1, plt2
    end

    return Plots.plot(plt1, plt2, layout=2, size=(2 * (xe-xs), ye-ys) .* size_mult)
end

rectangle((xs, xe), (ys, ye)) = Plots.Shape([xs, xe, xe, xs], [ys, ys, ye, ye])

function plot_comparison_for_cell(df_spatial::DataFrame, cell_id::Int, args...; cell1_col::Symbol=:cell, cell2_col::Symbol=:cell_paper, offset::Float64=0.1, kwargs...)
    sample_df = df_spatial[df_spatial[!, cell1_col] .== cell_id,:];
    paper_ids = setdiff(unique(sample_df[!, cell2_col]), [0]);

    if !isempty(paper_ids)
        sample_df = df_spatial[(df_spatial[!, cell1_col] .== cell_id) .| in.(df_spatial[!, cell2_col], Ref(paper_ids)),:];
    end;

    xc, yc = median(sample_df.x[sample_df[!, cell1_col] .== cell_id]), median(sample_df.y[sample_df[!, cell1_col] .== cell_id])
    xls, yls = [round.(Int, ((1. + offset) * s - offset * e, (1. + offset) * e - offset * s)) for (s,e) in [val_range(sample_df.x), val_range(sample_df.y)]];
    return plot_comparison_for_cell(df_spatial, xls, yls, args...; xc=xc, yc=yc, kwargs...)
end

function plot_comparison_for_cell(df_spatial::DataFrame, xls::Tuple{T, T}, yls::Tuple{T, T}, seg_arr::Matrix{<:Integer}, dapi_arr::Matrix{<:Real};
        size_mult::Float64=1.0, grid_alpha::Float64=0.0, ms::Float64=2.0, title="", center_mult::Float64=3.0, noise::Bool=false,
        xc::Union{Float64, Nothing}=nothing, yc::Union{Float64, Nothing}=nothing, kwargs...) where T <: Real

    xls, yls = max.(xls, 1), max.(yls, 1)
    xls = min.(xls, size(dapi_arr, 2))
    yls = min.(yls, size(dapi_arr, 1))

    seg_labels = seg_arr[yls[1]:yls[2], xls[1]:xls[2]];

    plts = plot_subset(df_spatial, dapi_arr, xls, yls; size_mult=size_mult, build_panel=false, grid_alpha=grid_alpha, ms=ms, noise=noise, kwargs...);

    paper_polys = [Plots.Shape(pg[:, 1], pg[:, 2]) for pg in extract_polygons_from_label_grid(copy(seg_labels'))]
    for plt in plts
        Plots.plot!(plt, paper_polys, fill=(0, 0.0), linewidth=1.5, alpha=0.75, linecolor="darkred", legend=:none);
        if xc !== nothing
            Plots.scatter!([xc - xls[1]], [yc - yls[1]], color="black", ms=center_mult*ms)
        end
    end;

    Plots.plot(plts..., layout=2, size=(plts[1].attr[:size] .* (2, 1)), title=title)
end

function joint_ordering(matrices::Matrix{T}...) where T <: Real
    joint_mtx = hcat(matrices...)
    cell_dists = Symmetric(1 .- cor(joint_mtx));
    cell_ord = Clustering.hclust(cell_dists, linkage=:ward).order;

    cell_ords = [Int[] for m in matrices]
    cum_sum = 0
    for i in 1:length(matrices)
        cell_ords[i] = cell_ord[(cell_ord .<= (cum_sum + size(matrices[i], 2))) .& (cell_ord .> cum_sum)] .- cum_sum
        cum_sum += size(matrices[i], 2)
    end

    gene_dists = Symmetric(1 .- cor(joint_mtx'));
    gene_ord = Clustering.hclust(gene_dists, linkage=:ward).order;

    return cell_ords, gene_ord
end