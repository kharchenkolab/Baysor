import PlotlyJS
using StatsBase
PLY = PlotlyJS

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

function plot_cur_state(bm_data::BmmData, subs_mask::BitArray{1}, clust_per_mol, cell_ids::Union{Vector{Int}, Nothing}=nothing;
                        gene_names=nothing, n_sigdigits::Int = 2, colorby="clust", line_width_mult::Float64 = 1e-4)
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
    end

    # Edges
    t_pos_data = position_data(bm_data.x);
    line_traces = vcat([[PLY.scatter(x=t_pos_data[1, [p1, p2]], y=t_pos_data[2, [p1, p2]], mode="lines", marker_color="black",
                    showlegend=false, hoverinfo="none", line_width=w * line_width_mult)
                for (w, p2) in zip(bm_data.adjacent_weights[p1], bm_data.adjacent_points[p1])]
            for p1 in mol_inds]...);

    # Text
    clust_per_mol = clust_per_mol[subs_mask]

    dens_comp = [round.(pdf.(Ref(bm_data.components[i].composition_params), c_subs.gene), sigdigits=n_sigdigits) for i in cell_ids];
    dens_pos = [round.(pdf.(Ref(bm_data.components[i].position_params), c_subs.x, c_subs.y), sigdigits=n_sigdigits) for i in cell_ids];

    if gene_names !== nothing
        c_subs.gene = gene_names[c_subs.gene]
    end

    plot_text = String[]
    for (i, (g, cl, cell)) in enumerate(zip(c_subs.gene, clust_per_mol, c_subs.cell))
        text = "Molecule: $i; Gene: $g;<br>Cluster: $cl; Cell: $cell;<br>Position:<br>"
        for (ci, d) in zip(cell_ids, dens_pos)
            text *= "$ci: $(d[i]); "
        end

        text *= "<br>Composition:<br>"
        for (ci, d) in zip(cell_ids, dens_comp)
            text *= "$ci: $(d[i]); "
        end

        adj_classes, adj_weights = adjacent_component_weights(bm_data.assignment, bm_data.adjacent_points[mol_inds[i]], bm_data.adjacent_weights[mol_inds[i]])[1:2]
        s_ids = sortperm(adj_classes)
        adj_classes, adj_weights = adj_classes[s_ids], adj_weights[s_ids]

        text *= "<br>Weights:<br>"
        for (ci, d) in zip(adj_classes, adj_weights)
            text *= "$ci: $(round(d, digits=1)); "
        end

        denses = adj_weights .* [c.prior_probability * pdf(c, c_subs.x[i], c_subs.y[i], c_subs.gene[i]) for c in bm_data.components[adj_classes]]
        denses ./= sum(denses)

        text *= "<br>Probs:<br>"
        for (ci, d) in zip(adj_classes, denses)
            text *= "$ci: $(round(d, sigdigits=n_sigdigits)); "
        end

        push!(plot_text, text)
    end

    # Molecules
    marker_color = (colorby == "clust") ? denserank(clust_per_mol) : denserank(c_subs.cell)
    marker_symbol = (colorby == "clust") ? denserank(c_subs.cell) : denserank(clust_per_mol)
    mol_trace = PLY.scatter(x=c_subs.x, y=c_subs.y, marker=PLY.attr(size=5, color=marker_color, symbol=marker_symbol, colorscale="Jet"),
        text=plot_text, mode="markers")

    return PLY.plot(vcat(line_traces, mol_trace), PLY.Layout(width=750, height=750, hovermode="closest"))
end

function plot_sampling_dynamic(bm_data::BmmData, subs_mask::BitArray{1}, clust_per_mol::Vector;
                               n_steps::Int=length(bm_data.tracer["assignment_history"]), step::Int=1)
    # See https://plot.ly/python/reference/#layout-sliders
    plot_texts = [["Cluster: $cl;<br>Gene: $g;<br>Cell: $cell" for (g, cl, cell) in zip(bm_data.x.gene, clust_per_mol[subs_mask], assignment)]
        for assignment in bm_data.tracer["assignment_history"]]

    steps = 1:step:n_steps
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
        hovermode="closest",
        width=700,
        height=700
    );

    return PLY.plot(PLY.scatter(x=bm_data.x.x, y=bm_data.x.y, text=plot_texts[1], mode=:markers), layout);
end

function set_to_step(bm_data::BmmData, iter::Int)
    bd_res = deepcopy(bm_data)
    bd_res.components = deepcopy(bd_res.tracer["component_history"][iter])

    # guid_to_luid = Dict(c.guid => i for (i, c) in enumerate(bd_res.components))
    # bd_res.assignment = get.(Ref(guid_to_luid), bd_res.tracer["assignment_history"][iter], -1)
    bd_res.assignment = deepcopy(bd_res.tracer["assignment_history"][iter])

    return bd_res
end
