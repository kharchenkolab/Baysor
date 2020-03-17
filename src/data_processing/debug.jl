import PlotlyJS
import PlotlyJS: PlotlyBase
using StatsBase
PLY = PlotlyJS

# TODO: Remove with PlotlyJS dependency?
function savehtml(plot, path::String)
    json = PLY.json(plot);
    html_content = """
    <html>
    <head>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <div id="4f98b041-213a-409c-9b89-fe52f87c5102" class="plotly-graph-div" style="float:left;margin-right:40px"></div>

        <script>
            window.PLOTLYENV=window.PLOTLYENV || {};
            window.PLOTLYENV.BASE_URL="https://plot.ly";
            Plotly.newPlot('4f98b041-213a-409c-9b89-fe52f87c5102', $json)
        </script>
    </body>
    </html>
    """;

    open(path, "w") do f
        print(f, html_content)
    end;
end

"""
    Requires running `bmm!` with `trace_components=true`
    Example:
        bds_it = Baysor.set_to_step(bmd_subs, 2);
        plot_cell_ids=[12, 13, 14, 15, 16]
        Baysor.plot_cur_state(bds_it, in.(bds_it.assignment, Ref(plot_cell_ids)), cell_ids=plot_cell_ids,
            clust_per_mol=bmd_subs.x.mol_clust, gene_names=gene_names, line_width_mult=1e-3)
"""
function plot_cur_state(bm_data::BmmData, subs_mask::BitArray{1}=trues(size(bm_data.x, 1)); cell_ids::Union{Vector{Int}, Nothing}=nothing,
                        clust_per_mol, adj_classes_global::Dict{Int, Vector{Int}}=Dict{Int, Vector{Int}}(), gene_names=nothing,
                        n_sigdigits::Int = 2, colorby="clust", line_width_mult::Float64 = 1e-4)
    c_subs = nothing
    mol_inds = nothing

    if length(subs_mask) == size(bm_data.x, 1)
        c_subs = deepcopy(bm_data.x[subs_mask,:])
        c_subs.cell = bm_data.assignment[subs_mask]
        mol_inds = findall(subs_mask)
    elseif sum(subs_mask) == size(bm_data.x, 1)
        c_subs = deepcopy(bm_data.x)
        c_subs.cell = bm_data.assignment
        mol_inds = 1:length(bm_data.assignment)
    else
        error("Can't fold sub_mask and bm_data.x together")
    end

    if cell_ids === nothing
        cell_ids = sort(unique(bm_data.assignment[mol_inds]))
        cell_ids = cell_ids[cell_ids .> 0]
    end

    # Edges
    t_pos_data = position_data(bm_data.x);

    line_traces::Array{PlotlyBase.GenericTrace, 1} = vcat([[PLY.scatter(x=t_pos_data[1, [p1, p2]], y=t_pos_data[2, [p1, p2]], mode="lines", marker_color="black",
                    showlegend=false, hoverinfo="none", line_width=w * line_width_mult)
                for (w, p2) in zip(bm_data.adjacent_weights[p1], bm_data.adjacent_points[p1])]
            for p1 in mol_inds]...);

    # Text
    clust_per_mol = clust_per_mol[subs_mask]

    # dens_comp = [round.(pdf.(Ref(bm_data.components[i].composition_params), c_subs.gene), sigdigits=n_sigdigits) for i in cell_ids];
    dens_comp = [round.(pdf_comp.(Ref(bm_data.components[i]), c_subs.gene), sigdigits=n_sigdigits) for i in cell_ids];
    dens_pos = [round.(pdf.(Ref(bm_data.components[i].position_params), c_subs.x, c_subs.y), sigdigits=n_sigdigits) for i in cell_ids];

    plot_text = String[]
    for (i, (g, cl, cell)) in enumerate(zip(c_subs.gene, clust_per_mol, c_subs.cell))
        gene_text =  (gene_names === nothing) ? "Gene ID: $g" : "Gene ID: $g, Gene: $(gene_names[g])"
        text = "Molecule: $i; $gene_text;<br>Cluster: $cl; Cell: $cell;<br>Position:<br>"
        for (ci, d) in zip(cell_ids, dens_pos)
            text *= "$ci: $(d[i]); "
        end

        text *= "<br>Composition:<br>"
        for (ci, d) in zip(cell_ids, dens_comp)
            text *= "$ci: $(d[i]); "
        end

        adj_classes = Int[]
        adj_weights = Float64[]
        adjacent_component_weights!(adj_weights, adj_classes, Dict{Int, Float64}(), bm_data.assignment, bm_data.adjacent_points[mol_inds[i]], bm_data.adjacent_weights[mol_inds[i]])
        adj_global = get(adj_classes_global, i, Int[]);
        if length(adj_global) > 0
            append!(adj_classes, adj_global)
            append!(adj_weights, ones(length(adj_global)) .* bm_data.real_edge_weight)
        end

        s_ids = sortperm(adj_classes)
        adj_classes, adj_weights = adj_classes[s_ids], adj_weights[s_ids]

        text *= "<br>Weights:<br>"
        for (ci, d) in zip(adj_classes, adj_weights)
            text *= "$ci: $(round(d, digits=1)); "
        end

        denses = adj_weights .* [c.prior_probability * pdf(c, c_subs.x[i], c_subs.y[i], c_subs.gene[i]) for c in bm_data.components[adj_classes]]
        denses ./= sum(denses)

        text *= "<br>Prior probabilities:<br>"
        for (ci, c) in zip(adj_classes, bm_data.components[adj_classes])
            text *= "$ci: $(round(c.prior_probability, sigdigits=n_sigdigits)); "
        end

        text *= "<br>Probs:<br>"
        for (ci, d) in zip(adj_classes, denses)
            text *= "$ci: $(round(d, sigdigits=n_sigdigits)); "
        end

        push!(plot_text, text)
    end

    # Molecules
    marker_color = (colorby == "clust") ? denserank(clust_per_mol) : denserank(c_subs.cell)
    marker_symbol = (colorby == "clust") ? denserank(c_subs.cell) : denserank(clust_per_mol)
    mol_trace::PlotlyBase.GenericTrace = PLY.scatter(x=c_subs.x, y=c_subs.y, marker=PLY.attr(size=5, color=marker_color, symbol=marker_symbol, colorscale="Jet"),
        text=plot_text, mode="markers")

    return PLY.plot(vcat(line_traces, mol_trace), PLY.Layout(width=750, height=750, hovermode="closest"))
end

function plot_sampling_dynamic(bm_data::BmmData, subs_mask::BitArray{1}, clust_per_mol::Vector;
                               start_step::Int=1, n_steps::Int=length(bm_data.tracer["assignment_history"]), step::Int=1)
    # See https://plot.ly/python/reference/#layout-sliders
    plot_texts = [["Cluster: $cl;<br>Gene: $g;<br>Cell: $cell" for (g, cl, cell) in zip(bm_data.x.gene, clust_per_mol[subs_mask], assignment)]
        for assignment in bm_data.tracer["assignment_history"]]

    steps = start_step:step:n_steps
    layout = PLY.Layout(
        sliders=[PLY.attr(
            steps=[PLY.attr(method="restyle",
                        args=[Dict(
                                "marker" => PLY.attr(size=5, color=denserank(cols), symbol=denserank(clust_per_mol[subs_mask]), colorscale="Jet"),
                                "text" => [text],
                                "visible" => vcat(true, steps .== i)
                            )], label=i)
                    for (i, cols, text) in zip(steps, bm_data.tracer["assignment_history"][steps], plot_texts[steps])],
            active=0,
            currentvalue_prefix="Step: "
        )],
        updatemenus=[PLY.attr(
            type="buttons",
            buttons=[PLY.attr(
                args=[nothing],
                label="Play",
                method="restyle"
            )]

        )],
        hovermode="closest",
        width=700,
        height=700
    );

    return PLY.plot(PLY.scatter(x=bm_data.x.x, y=bm_data.x.y, text=plot_texts[1], mode=:markers), layout);
end

"""
    Example:
        polygons_per_iter = Baysor.boundary_polygons.(Ref(bmd_subs.x), bmd_subs.tracer[:assignment_history], grid_step=0.5,
            min_molecules_per_cell=5, bandwidth=0.5);
        Baysor.plot_sampling_dynamic(bmd_subs, polygons_per_iter; step_ids=1:50, width=600, height=600, marker_color=bmd_subs.x.mol_clust)
"""
function plot_sampling_dynamic(bm_data::BmmData, polygons_per_iter::Array; step_ids::T where T <: AbstractArray{Int}=1:length(polygons_per_iter),
                               min_polygon_size::Int=5, width::Int=600, height::Int=600, marker_color)
    function plotly_polygon(polygon::Matrix; kwargs...)
        polygon = vcat(polygon, polygon[1,:]')
        return PLY.scatter(x=polygon[:,1], y=polygon[:,2], mode="lines", marker_color="black", showlegend=false, hoverinfo="none"; kwargs...)
    end

    polygon_shapes = [plotly_polygon.(filter(p -> size(p, 1) >= min_polygon_size, ps); visible=(i == 1))
        for (i, ps) in enumerate(polygons_per_iter[step_ids])];

    n_pols_total = sum(length.(polygon_shapes));
    n_pols_cum = 1
    is_visible = []

    for n_pols in length.(polygon_shapes)
        cur_mask = falses(n_pols_total)
        cur_mask[n_pols_cum:(n_pols_cum + n_pols - 1)] .= true
        n_pols_cum += n_pols
        push!(is_visible, cur_mask)
    end

    assignment_history = bm_data.tracer[:assignment_history][step_ids];
    texts = [["molecule: $mid, cell: $cid<br>x: $x, y: $y" for (mid, cid, x, y) in zip(1:size(bm_data.x, 1), ah, bm_data.x.x, bm_data.x.y)] for ah in assignment_history]

    mol_trace = PLY.scatter(x=bm_data.x.x, y=bm_data.x.y, text=texts[1], marker=PLY.attr(size=5, color=marker_color, colorscale="Jet"),
        showlegend=false, mode="markers")

    layout = PLY.Layout(
        sliders=[PLY.attr(
            # steps=vcat(PLY.attr(method="restyle", args=[Dict("visible" => vcat(true, falses(n_pols_total)))], label=0),
            #     [PLY.attr(method="restyle", args=[Dict("visible" => vcat(true, v))], label=i) for (i, v) in zip(step_ids, is_visible)]
            # ),
            steps=[PLY.attr(method="restyle", args=[Dict("visible" => vcat(true, v), "text" => [txt])], label=i) for (i, v, txt) in zip(step_ids, is_visible, texts)],
            active=0,
            currentvalue_prefix="Step: "
        )],
        hovermode="closest",
        width=width,
        height=height
    );

    return PLY.plot(vcat(mol_trace, polygon_shapes...), layout);
end

function set_to_step(bm_data::BmmData, iter::Int)
    bd_res = deepcopy(bm_data)
    bd_res.components = deepcopy(bd_res.tracer[:component_history][iter])

    # guid_to_luid = Dict(c.guid => i for (i, c) in enumerate(bd_res.components))
    # bd_res.assignment = get.(Ref(guid_to_luid), bd_res.tracer["assignment_history"][iter], -1)
    bd_res.assignment = deepcopy(bd_res.tracer[:assignment_history][iter])

    return bd_res
end

function plot_expression_vectors(vecs...; gene_names::Vector{String}, min_expr::Float64=0.05, alpha::Float64=0.5, fontsize::Int=5, text_offset::Float64=0.005, kwargs...)
    p = Plots.plot(;kwargs...)
    for v in vecs
        p = Plots.bar!(v, alpha=alpha)
    end

    y_vals = maximum(hcat(vecs...), dims=2) |> vec
    ann_genes = findall(y_vals .>= min_expr)
    p = Plots.annotate!(ann_genes, y_vals[ann_genes] .+ text_offset, Plots.text.(gene_names[ann_genes], fontsize))
    return p
end
