using DataFrames

struct InitialParams
    centers::Array{Float64, 2}
    stds::Array{Array{Float64,2},1};
    assignment::Array{Int64,1};
    n_comps::Int;

    function InitialParams(centers::Array{T, 2} where T <: Real, stds::Union{Array{Array{Float64,2},1}, Real}, assignment::Array{Int64,1})
        @assert size(centers, 1) == maximum(assignment)
        if isa(stds, Real)
            stds = [[stds 0.; 0. stds] for i in 1:size(centers, 1)]
        else
            @assert size(centers, 1) == size(stds, 1)
        end

        return new(centers, adjust_cov_matrix.(stds), assignment, size(centers, 1))
    end
end

mutable struct BmmData
    components::Array{Component, 1};

    x::DataFrame;
    position_data::Array{Float64, 2};
    composition_data::Array{Int, 1};

    adjacent_points::Array{Array{Int, 1}, 1};
    distribution_sampler::Component;
    assignment::Array{Int, 1};
    position_knn_tree::KDTree;
    knn_neighbors::Array{Array{Int, 1}, 1};

    gene_probs_given_single_transcript::Array{Float64, 2};

    tracer::Dict{String, Any};

    function BmmData(components::Array{Component, 1}, x::DataFrame, adjacent_points::Array{Array{Int, 1}, 1}, distribution_sampler::Component,
                     assignment::Array{Int, 1}; k_neighbors::Int=20)
        for c in components
            c.n_samples = 0
        end

        p_data = position_data(x)
        position_knn_tree = KDTree(p_data)
        knn_neighbors = knn(position_knn_tree, p_data, k_neighbors)[1]

        n_genes = maximum(composition_data(x))

        self = new(components, deepcopy(x), p_data, composition_data(x), adjacent_points, deepcopy(distribution_sampler), zeros(Int, length(assignment)),
                   position_knn_tree, knn_neighbors, ones(n_genes, n_genes) ./ n_genes, Dict{String, Any}())

        for (i, l) in enumerate(assignment)
            assign!(self, i, l)
        end

        return self
    end
end

num_of_molecules_per_cell(data::BmmData) = count_array(data.assignment .+ 1, max_value=length(data.components) + 1)[2:end]

function merge_bm_data(bmm_data_arr::Array{BmmData, 1})
    @assert length(bmm_data_arr) > 0

    x = vcat([bd.x for bd in bmm_data_arr]...)
    components = vcat([bd.components for bd in bmm_data_arr]...)

    adjacent_points = Array{Int64,1}[]
    ap_offset = 0;
    for bd in bmm_data_arr
        append!(adjacent_points, [ap .+ ap_offset for ap in bd.adjacent_points])
        ap_offset += size(bd.x, 1)
    end

    assignments = Array{Int64,1}[]
    assignment_offset = 0;
    for bd in bmm_data_arr
        cur_assignment = bd.assignment .+ assignment_offset
        cur_assignment[bd.assignment .== 0] .= 0
        push!(assignments, cur_assignment)
        assignment_offset += length(bd.components)
    end

    k_neighbors=length(bmm_data_arr[1].knn_neighbors[1])

    res = BmmData(components, x, adjacent_points, bmm_data_arr[1].distribution_sampler, vcat(assignments...); k_neighbors=k_neighbors)
    res.tracer = merge_tracers([bd.tracer for bd in bmm_data_arr])

    return res
end