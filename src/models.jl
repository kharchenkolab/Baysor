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
                   position_knn_tree, knn_neighbors, ones(n_genes, n_genes) ./ n_genes)

        for (i, l) in enumerate(assignment)
            assign!(self, i, l)
        end

        return self
    end
end

num_of_molecules_per_cell(data::BmmData) = count_array(data.assignment .+ 1, max_value=length(data.components) + 1)[2:end]