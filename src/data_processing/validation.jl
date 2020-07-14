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
        "n_cells" => [sum(count_array(x, drop_zero=true) .>= min_molecules_per_cell) for x in values(assignments)]
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

function plot_subset(df_spatial::DataFrame, dapi_arr::Union{Matrix{<:Real}, Nothing}, (xs, xe), (ys, ye); polygons::Union{Bool, Vector{Matrix{Float64}}}=true, ms=2.0, alpha=0.2, min_molecules_per_cell::Int=1,
        grid_step::Float64=5.0, bandwidth::Float64=grid_step, cell_col::Symbol=:cell, dapi_alpha=0.9, polygon_line_width::T1 where T1 <: Real=2, polygon_alpha::Float64=0.4, dens_threshold::Float64=1e-5,
        noise::Bool=true, size_mult=1/3, plot_raw_dapi::Bool=true, color_col::Symbol=:color, annotation_col::Union{Symbol, Nothing}=nothing, build_panel::Bool=true, grid_alpha::Float64=0.5, ticks=false,
        grid::Bool=true, kwargs...)
    df_subs = @where(df_spatial, :x .>= xs, :x .<= xe, :y .>= ys, :y .<= ye);

    if (typeof(polygons) == Bool)
        if polygons
            polygons = boundary_polygons(df_subs, df_subs[!, cell_col], grid_step=grid_step, min_molecules_per_cell=min_molecules_per_cell, bandwidth=bandwidth, dens_threshold=dens_threshold)
        else
            polygons = Matrix{Float64}[]
        end
    end

    # xs, xe, ys, ye = round.(Int, [minimum(df_subs.x), maximum(df_subs.x), minimum(df_subs.y), maximum(df_subs.y)]);
    xs, xe, ys, ye = round.(Int, [xs, xe, ys, ye]);

    xticks_vals = range(0, xe-xs, length=5)
    yticks_vals = range(0, ye-ys, length=5)

    yticks = xticks = Float64[]
    if ticks
        xticks = (xticks_vals, ["$f" for f in round.(Int, range(xs, xe, length=5))])
        yticks = (yticks_vals, ["$f" for f in round.(Int, range(ys, ye, length=5))])
    end

    plot_size = ((xe-xs), ye-ys) .* size_mult
    is_noise = noise ? (df_subs[!, cell_col] .== 0) : nothing
    annotation = (annotation_col === nothing) ? nothing : df_subs[!, annotation_col]

    if dapi_arr === nothing
        return plot_cell_borders_polygons(df_subs, polygons; color=df_subs[!, color_col], ms=ms, alpha=alpha, offset=(-xs, -ys), size=plot_size,
            polygon_line_width=polygon_line_width, polygon_alpha=polygon_alpha, is_noise=is_noise, noise_kwargs=Dict(:ms => 1.0), annotation=annotation, kwargs...)
    end

    dapi_subs = dapi_arr[ys:ye, xs:xe]
    plt1 = Plots.heatmap(dapi_subs, color=:grayC, colorbar=:none, size=plot_size,
                alpha=dapi_alpha, format=:png, legend=:none, xticks=xticks, yticks=yticks)

    plot_cell_borders_polygons!(df_subs, polygons; color=df_subs[!, color_col], ms=ms, alpha=alpha, offset=(-xs, -ys),
        polygon_line_width=polygon_line_width, polygon_alpha=polygon_alpha, is_noise=is_noise, noise_kwargs=Dict(:ms => 1.0), annotation=annotation, kwargs...)

    if grid
        Plots.vline!(xticks_vals, color="black", alpha=grid_alpha, label="")
        Plots.hline!(yticks_vals, color="black", alpha=grid_alpha, label="")
    end

    if !plot_raw_dapi
        return plt1
    end

    # Possible colorschemes: tarn, diff, lime_grad, thermal
    plt2 = Plots.heatmap(dapi_subs, color=:delta, colorbar=:none, alpha=0.9, format=:png, xticks=xticks, yticks=yticks, legend=:none, size=plot_size)

    Plots.plot!([Plots.Shape(pg[:,1] .- xs, pg[:,2] .- ys) for pg in polygons],
        fill=(0, 0.0), linewidth=polygon_line_width, linecolor="black", alpha=polygon_alpha, label="", xlims=(0, (xe-xs)), ylims=(0, (ye-ys)));

    if grid
        Plots.vline!(xticks_vals, color="black", alpha=grid_alpha)
        Plots.hline!(yticks_vals, color="black", alpha=grid_alpha)
    end

    if !build_panel
        return plt1, plt2
    end

    return Plots.plot(plt1, plt2, layout=2, size=(2 * (xe-xs), ye-ys) .* size_mult)
end

rectangle((xs, xe), (ys, ye)) = Plots.Shape([xs, xe, xe, xs], [ys, ys, ye, ye])

function cell_coord_frame(df_spatial::DataFrame, cell_id::Int; cell_col::Symbol=:cell, offset::Float64=0.1)
    sample_df = df_spatial[df_spatial[!, cell_col] .== cell_id,:];
    return [round.(Int, ((1. + offset) * s - offset * e, (1. + offset) * e - offset * s)) for (s,e) in [val_range(sample_df.x), val_range(sample_df.y)]];
end

function plot_comparison_for_cell(df_spatial::DataFrame, cell_id::Int, args...; cell1_col::Symbol=:cell, cell2_col::Symbol=:cell_dapi, offset::Float64=0.1, kwargs...)
    sample_df = df_spatial[df_spatial[!, cell1_col] .== cell_id,:];
    paper_ids = setdiff(unique(sample_df[!, cell2_col]), [0]);

    if !isempty(paper_ids)
        sample_df = df_spatial[(df_spatial[!, cell1_col] .== cell_id) .| in.(df_spatial[!, cell2_col], Ref(paper_ids)),:];
    end;

    xc, yc = median(sample_df.x[sample_df[!, cell1_col] .== cell_id]), median(sample_df.y[sample_df[!, cell1_col] .== cell_id])
    xls, yls = [round.(Int, ((1. + offset) * s - offset * e, (1. + offset) * e - offset * s)) for (s,e) in [val_range(sample_df.x), val_range(sample_df.y)]];
    return plot_comparison_for_cell(df_spatial, xls, yls, args...; xc=xc, yc=yc, cell_col=cell1_col, kwargs...)
end

function plot_comparison_for_cell(df_spatial::DataFrame, xls::Tuple{T, T}, yls::Tuple{T, T}, seg_arr::Union{Matrix{<:Integer}, Nothing},
        dapi_arr::Union{Matrix{<:Real}, Nothing}; paper_polys::Array{Matrix{Float64}, 1}=Matrix{Float64}[], polygon_line_width::Float64=2.0, polygon_alpha::Float64=0.5,
        size_mult::Float64=1.0, grid_alpha::Float64=0.0, ms::Float64=2.0, title="", center_mult::Float64=3.0, noise::Bool=false, paper_poly_color="darkred", paper_line_mult::Float64=1.0,
        xc::Union{Float64, Nothing}=nothing, yc::Union{Float64, Nothing}=nothing, plot_raw_dapi::Bool=true, kwargs...) where T <: Real

    xls, yls = max.(xls, 1), max.(yls, 1)

    if dapi_arr !== nothing
        xls = min.(xls, size(dapi_arr, 2))
        yls = min.(yls, size(dapi_arr, 1))
    end

    if seg_arr !== nothing
        paper_polys = extract_polygons_from_label_grid(copy(seg_arr[yls[1]:yls[2], xls[1]:xls[2]]'))
        paper_polys = [Plots.Shape(pg[:,1], pg[:,2]) for pg in paper_polys]
    else
        paper_polys = [Plots.Shape(pg[:,1] .- xls[1], pg[:,2] .- yls[1]) for pg in paper_polys]
    end

    if dapi_arr === nothing
        plot_raw_dapi = false
    end

    plts = plot_subset(df_spatial, dapi_arr, xls, yls; size_mult=size_mult, build_panel=false, grid_alpha=grid_alpha, ms=ms, noise=noise,
        polygon_line_width=polygon_line_width, polygon_alpha=polygon_alpha, plot_raw_dapi=plot_raw_dapi, kwargs...);
    if !plot_raw_dapi
        plts = [plts]
    end

    for plt in plts
        Plots.plot!(plt, paper_polys, fill=(0, 0.0), linewidth=(polygon_line_width*paper_line_mult), alpha=polygon_alpha, linecolor=paper_poly_color, legend=:none);
        if xc !== nothing
            Plots.scatter!([xc - xls[1]], [yc - yls[1]], color="black", ms=center_mult*ms)
        end
    end;

    if !plot_raw_dapi
        return plts[1]
    end

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

## Comparison of segmentations

function prepare_qc_df(df_spatial::DataFrame, cell_col::Symbol=:cell; min_area::T where T<:Real, min_molecules_per_cell::T2 where T2 <: Real,
        max_elongation::T3 where T3 <: Real = 30.0, dapi_arr::Union{Matrix{<:Real}, Nothing}=nothing)
    qc_per_cell = get_cell_qc_df(df_spatial, Int.(df_spatial[!, cell_col]), dapi_arr=dapi_arr);
    qc_per_cell[!, :cell_id] = 1:size(qc_per_cell, 1)
    qc_per_cell[!, :sqr_area] = sqrt.(qc_per_cell.area)
    return @where(qc_per_cell, :n_transcripts .>= min_molecules_per_cell, :sqr_area .>= min_area, :elongation .<= max_elongation)
end

hist_bins(vals::Vector{<:Real}...; n_bins::Int=100, min_val::T where T<: Real=0.0, m_quantile::T2 where T2<:Real=0.99, max_val_mult::Float64=1.0) =
    range(min_val, maximum([quantile(v, m_quantile) / m_quantile for v in vals if length(v) > 0])*max_val_mult, length=n_bins)

function plot_qc_comparison(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}; labels::Vector{String}=["Baysor", "DAPI"],
        max_quants::Vector{<:Real}=[0.9999, 0.99, 0.99, 0.999, 0.999], n_bins::Vector{<:Real}=[75, 100, 100, 100, 100], kwargs...)
    m_names = [:n_transcripts, :density, :elongation, :sqr_area, :mean_dapi_brightness]
    plot_titles = ["Num. of transcripts", "Density", "Elongation", "sqrt(Area)", "Mean DAPI brightness"]
    plots = Plots.Plot[]
    for (cs, xlab, mq, nb) in zip(m_names, plot_titles, max_quants, n_bins)
        if any(!(cs in names(qdf)) for qdf in qc_per_cell_dfs)
            continue
        end
        t_bins = hist_bins([qdf[!, cs] for qdf in qc_per_cell_dfs]..., m_quantile=mq, n_bins=nb)
        plt = Plots.plot(widen=false, xlabel=xlab, ylabel="Num. of cells", xlims=Baysor.val_range(t_bins))
        for (qdf, lab) in zip(qc_per_cell_dfs, labels)
            if length(unique(qdf[!, cs])) < 2
                continue
            end
            Plots.histogram!(qdf[!, cs], label=lab, bins=t_bins, alpha=0.6, kwargs...)
        end
        push!(plots, plt)
    end

    return Plots.plot(plots..., size=(400 * ceil(Int, length(plots) / 2), 600))
end

function match_assignments(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, df_spatial::DataFrame, cell_cols::Vector{Symbol}; kwargs...)
    assignments = [denserank(ifelse.(in.(df_spatial[!,s], Ref(Set(qdf.cell_id))), df_spatial[!,s], 0)) .- 1 for (s, qdf) in zip(cell_cols, qc_per_cell_dfs)];
    return match_assignments(assignments...; kwargs...)
end

function match_assignments(assignment1::Vector{<:Integer}, assignment2::Vector{<:Integer}; bin_match_frac::Float64=0.05)
    # Basic statistics
    contingency_table = sparse_counts(assignment1 .+ 1, assignment2 .+ 1);
    ci,ri,v = SparseArrays.findnz(contingency_table)
    m_norm1 = sparse(ci, ri, v ./ sum(contingency_table, dims=1)[ri])
    m_norm2 = sparse(ci, ri, v ./ sum(contingency_table, dims=2)[ci])

    noise_fracs = Vector.((m_norm2[2:end,1], m_norm1[1,2:end]));
    max_overlaps = Vector.((maximum(m_norm2[2:end,:], dims=2)[:], maximum(m_norm1[:,2:end], dims=1)[:]));

    # @time bin_match = (contingency_table[2:end, 2:end] ./ sum(contingency_table[2:end, 2:end], dims=1) .> bin_match_frac);
    bin_match = dropzeros!(sparse(ci, ri, m_norm1.nzval .> bin_match_frac))[2:end, 2:end];
    n_overlaps = (sum(bin_match, dims=2)[:], sum(bin_match, dims=1)[:])

    ctn_min = sparse(ci, ri, min.(m_norm1.nzval, m_norm2.nzval))
    match_fracs = range(0.0, 1.0, length=30)
    match_nums = [sum(any(ctn_min .>= mt, dims=1)) for mt in match_fracs];

    # Derivative statistics
    match_noise = [nf .> 0.5 for nf in noise_fracs]
    multiple_overlap = [((mo .< 0.7) .& (nf .< 0.25)) for (mo, nf) in zip(max_overlaps, noise_fracs)];
    @assert all([!any(mn .& mo) for (mn, mo) in zip(match_noise, multiple_overlap)])

    return (contingency=contingency_table, noise_fracs=noise_fracs, max_overlaps=max_overlaps, n_overlaps=n_overlaps, match_fracs=match_fracs,
        match_nums=match_nums, match_noise=match_noise, multiple_overlap=multiple_overlap)
end

function plot_matching_comparison(match_results::NamedTuple; labels::Vector{String}=["Baysor", "DAPI"], n_bins::Vector{<:Real}=[50, 30, 30], kwargs...)
    m_names = [:noise_fracs, :max_overlaps]
    plot_titles = ["Fraction of molecules, not assigned in other segmentation", "Fraction of molecules, matching to a single cell"]
    plots = Plots.Plot[]
    for (n, xlab, nb) in zip(m_names, plot_titles, n_bins)
        vals = match_results[n]
        t_bins = hist_bins(vals..., m_quantile=1.0, n_bins=nb, max_val_mult=1.02)
        plt = Plots.histogram(vals[1], label=labels[1], bins=t_bins, xlims=val_range(t_bins), widen=false,
            xlabel=xlab, ylabel="Num. of cells", kwargs...)
        Plots.histogram!(vals[2], label=labels[2], bins=t_bins, alpha=0.6)
        push!(plots, plt)
    end

    max_bin = maximum(maximum.(match_results.n_overlaps))
    plt = Plots.histogram(match_results.n_overlaps[1], label=labels[1], bins=0:0.5:max_bin, xticks=0:max_bin, widen=false,
        xlabel="Number of overlapping cells", ylabel="Num. of cells", kwargs...)
    Plots.histogram!(match_results.n_overlaps[2], label=labels[2], bins=0:0.5:max_bin, alpha=0.6)
    push!(plots, plt)

    plt = Plots.plot(match_results.match_fracs, match_results.match_nums, xlabel="Minimal fraction of matching molecules", ylabel="Num. of matching cells",
        xlims=(0.0, 1.0), ylims=(0, maximum(match_results.match_nums) * 1.05), widen=false, legend=:none, lw=2.0)
    push!(plots, plt)

    return Plots.plot(plots..., layout=(2, 2), size=(900, 600))
end

function build_statistics_df(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, match_res::NamedTuple, df_spatial::DataFrame; labels::Vector{Symbol}=[:Baysor, :DAPI], sigdigits::Int=3)
    stat_df = DataFrame(Dict(:Metric => ["Num. cells", "Total num. molecules", "Fraction of noise", "Fraction of matching cells", "Fraction of matching to noise", "Fraction of multiple overlap"]))
    for i in 1:length(qc_per_cell_dfs)
        stat_df[!, labels[i]] = Any[
            size(qc_per_cell_dfs[i], 1), sum(qc_per_cell_dfs[i].n_transcripts), round(1 .- sum(qc_per_cell_dfs[i].n_transcripts) / size(df_spatial, 1), sigdigits=sigdigits),
            round(mean(.!match_res.match_noise[i] .& (.!match_res.multiple_overlap[i])), sigdigits=sigdigits), round(mean(match_res.match_noise[i]), sigdigits=sigdigits),
            round(mean(match_res.multiple_overlap[i]), sigdigits=sigdigits)
        ]
    end

    return stat_df
end

estimate_embedding(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, cell_cols::Vector{Symbol}; kwargs...) =
    estimate_embedding([convert_segmentation_to_counts(df_spatial.gene, df_spatial[!, cq])[:, qdf.cell_id] for (cq, qdf) in zip(cell_cols, qc_per_cell_dfs)]...; kwargs...)

function estimate_embedding(count_matrices::Matrix{<:Real}...; joint::Bool=true, kwargs...)
    cm_merged = hcat(count_matrices...);
    cm_merged = cm_merged ./ sum(cm_merged, dims=1)
    ids_per_mat = split_ids(vcat([repeat([i], inner=size(cm, 2)) for (i,cm) in enumerate(count_matrices)]...))
    umap_merged = fit(UmapFit, cm_merged; kwargs...);
    return [umap_merged.embedding[:, ids] for ids in ids_per_mat], umap_merged, cm_merged
end

function plot_qc_embeddings(qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, match_res::NamedTuple, embeddings::Array{Matrix{Float64}, 1}; log_colors::Bool=true, size=(1000, 1000),
        labels::Vector{String}=["Baysor", "DAPI"], legend::Symbol=:best, ms::Float64=3.0, alpha1::Float64=0.25, plot_num_transcripts::Bool=true, layout=nothing,
        plot_partial::Bool=false, palette=Plots.distinguishable_colors(4, cchoices=[80, 100]), build_figure::Bool=true, kwargs...)
    plts = Plots.Plot[]
    c_vals = [qdf.n_transcripts for qdf in qc_per_cell_dfs]
    if log_colors
        c_vals = [log.(cv) for cv in c_vals]
    end

    factor_per_cell = String[]

    if plot_partial
        factor_per_cell = [ifelse.(match_res.multiple_overlap[i], "Multiple", ifelse.(match_res.match_noise[i], "No",
            ifelse.(match_res.max_overlaps[i] .< 0.7, "Partial", "Full"))) for i in 1:length(qc_per_cell_dfs)];
    else
        factor_per_cell = [ifelse.(match_res.multiple_overlap[i], "Multiple", ifelse.(match_res.match_noise[i], "No", "Yes")) for i in 1:length(qc_per_cell_dfs)];
    end

    for i in 1:length(qc_per_cell_dfs)
        plt = Plots.plot(format=:png, size=(600, 600), legend=legend, bg_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), title="$(labels[i]), matching", legend_title="Match")
        for (ci,ids) in enumerate(split_ids(denserank(factor_per_cell[i])))
            Plots.scatter!(embeddings[i][1,ids], embeddings[i][2,ids], label=factor_per_cell[i][ids[1]], ms=ms, alpha=alpha1, markerstrokewidth=0; color=palette[ci], kwargs...)
        end
        push!(plts, plt)

        plot_num_transcripts || continue

        colors = map_to_colors(c_vals[i], lims=val_range(vcat(c_vals...)))[:colors]
        plt = Plots.scatter(embeddings[i][1,:], embeddings[i][2,:], color=colors, ms=ms, markerstrokewidth=0, format=:png, size=(600, 600), legend=false,
            title="$(labels[i]), num. transcripts")
        push!(plts, plt)
    end;

    if !build_figure
        return plts
    end

    return Plots.plot(plts..., size=size, format=:png, layout=something(layout, length(plts)))
end

plot_expression_vec_comparison(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}}, cell_cols::Vector{Symbol}, gene_names; labels=["Baysor", "DAPI"], kwargs...) =
    plot_expression_vectors([count_array(vcat(split(df_spatial.gene, df_spatial[!, cs], drop_zero=true)[qdf.cell_id]...)) for (cs, qdf) in zip(cell_cols, qc_per_cell_dfs)]...; gene_names=gene_names, labels=labels, kwargs...)

function estimate_non_matching_part_correlation(df_spatial::DataFrame, qc_per_cell_dfs::Union{Array{DataFrame, 1}, Tuple{DataFrame, DataFrame}},
        match_res::NamedTuple; cell_cols::Vector{Symbol}=[:cell, :cell_dapi], rev::Bool=false, max_overlap_borders::Tuple{Float64, Float64}=(0.25, 0.75))
    data_id = rev ? 2 : 1
    if rev
        cell_cols = reverse(cell_cols)
    end
    part_overlap_mask = (match_res.max_overlaps[data_id] .> max_overlap_borders[1]) .& (match_res.max_overlaps[data_id] .< max_overlap_borders[2]);
    cur_overlaps = match_res.max_overlaps[data_id][part_overlap_mask]
    partially_matching_ids = qc_per_cell_dfs[data_id].cell_id[part_overlap_mask];
    match_cors = Float64[]

    n_genes = maximum(df_spatial.gene)
    for t_cell_id in partially_matching_ids
        t_df = df_spatial[df_spatial[!, cell_cols[1]] .== t_cell_id,:];
        t_mcp = mode(t_df[!, cell_cols[2]][t_df[!, cell_cols[2]] .!= 0]);

        match_expr = prob_array(t_df.gene[t_df[!, cell_cols[2]] .== t_mcp], max_value=n_genes);
        non_match_expr = prob_array(t_df.gene[t_df[!, cell_cols[2]] .!= t_mcp], max_value=n_genes);

        push!(match_cors, cor(non_match_expr, match_expr))
    end

    return match_cors, partially_matching_ids, cur_overlaps, part_overlap_mask
end

function correlation_effect_size_df(datasets::NamedTuple, part_cor_key::Symbol, cell_col_labels::Vector{String}, cell_col_ids::Vector{Int})
    dfs = DataFrame[]
    for i in 1:2
        for d in datasets
            (part_cor_key in keys(d)) || continue

            corr_res = d[part_cor_key][i]
            qc_df = d[:qc_per_cell_dfs][cell_col_ids[i]]
            corr = corr_res[1];
            ord = sortperm(corr);
            corr = corr[ord]

            mismatch_frac = 1 .- corr_res[3][ord];
            n_trans_per_cell = qc_df.n_transcripts[corr_res[4]][ord];

            push!(dfs, DataFrame(
                :Correlation => corr,
                :MolFrac => cumsum(mismatch_frac .* n_trans_per_cell) ./ sum(qc_df.n_transcripts),
                :CellFrac => (1:length(corr)) ./ size(qc_df, 1),
                :Dataset => d[:name],
                :Segmentation => cell_col_labels[cell_col_ids[i]]
            ))
        end
    end
    return vcat(dfs...)
end
