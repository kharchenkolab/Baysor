using StatsBase: countmap

function get_best_molecule_assignment(adjacent_points::Array{Int,1}, adj_freqs::Dict{Int, Float64}, assignment::Array{Int, 1})::Int
    sum_freq_per_clust = Dict{Int, Float64}();

    for i in adjacent_points
        cur_assignment = assignment[i]
        if cur_assignment != 0
            sum_freq_per_clust[cur_assignment] = get(sum_freq_per_clust, cur_assignment, 0) + get(adj_freqs, i, 0)
        end
    end

    if length(sum_freq_per_clust) == 0
        return 0
    end

    return collect(keys(sum_freq_per_clust))[findmax(collect(values(sum_freq_per_clust)))[2]]
end

function extract_neighborhood_from_assignment(assignment::Array{Int, 1})
    neighbors_by_molecule = [Int[] for i in 1:length(assignment)]
    for ids in split(collect(1:length(assignment)), assignment .+ 1)[2:end]
        for i in ids
            neighbors_by_molecule[i] = copy(ids)
        end
    end

    return neighbors_by_molecule
end

function estimate_adjacent_point_frequencies(assignment_history::Array{Array{Array{Int,1},1},1}, n_stable_iters::Int)
    adj_freqs = [countmap(vcat([h[i] for h in assignment_history[end-n_stable_iters+1:end]]...)) for i in 1:length(assignment_history[1])];
    return [Dict(k => v / n_stable_iters for (k, v) in freqs) for freqs in adj_freqs]
end

function estimate_uncertainty(neighbors_by_molecule::Array{Array{Int, 1}}, adj_freqs::Array{Dict{Int, Float64}})
    confidence = [(length(neighb_ids) > 0 ? mean(Float64[get(dict, i, 0) for i in neighb_ids]) : 0) for (dict, neighb_ids) in zip(adj_freqs, neighbors_by_molecule)];
    return 1 .- confidence
end

estimate_uncertainty(neighbors_by_molecule::Array{Array{Int, 1}}, assignment_history::Array{Array{Array{Int,1},1},1}, n_stable_iters::Int) =
    estimate_uncertainty(neighbors_by_molecule, estimate_adjacent_point_frequencies(assignment_history, n_stable_iters))

function reassign_molecules_with_history(bm_data, assignment_history::Array{Array{Array{Int,1},1},1}, n_stable_iters::Int;
                                         min_molecules_per_cell::Int = 10, max_iters::Int=100, outlier_confidence_threshold::Float64=0.25,
                                         return_uncertainty::Bool=true)
    n_stable_iters = min(n_stable_iters, length(assignment_history))
    adj_freqs = estimate_adjacent_point_frequencies(assignment_history, n_stable_iters)
    freq_per_current_neighbors = [[dict[i] for i in hist] for (dict, hist) in zip(adj_freqs, assignment_history[end])];

    assignment::Array{Int, 1} = copy(bm_data.assignment);
    assignment[length.(freq_per_current_neighbors) .< min_molecules_per_cell] .= 0;

    for iter in 1:max_iters
        prev_assignment = copy(assignment)

        for mol_id in sortperm(mean.(freq_per_current_neighbors))
            assignment[mol_id] = get_best_molecule_assignment(bm_data.adj_list.ids[mol_id], adj_freqs[mol_id], assignment)
        end

        if all(prev_assignment .== assignment)
            break
        end
    end

    is_outlier = [sum(collect(values(d)) .> outlier_confidence_threshold) for d in adj_freqs] .< min_molecules_per_cell;
    assignment[is_outlier] .= 0;

    if !return_uncertainty
        return assignment
    end

    uncertainty = estimate_uncertainty(extract_neighborhood_from_assignment(assignment), adj_freqs)

    return assignment, uncertainty
end
