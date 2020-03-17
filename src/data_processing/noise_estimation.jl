function find_noise_border(dists_to_top::Vector{Float64}, dists_to_bottom::Vector{Float64}, mol_order::Vector{Int};
        is_reverse::Bool=false, passed_frac_threshold::Float64=0.95)::Int
    dists1, dists2 = dists_to_top, dists_to_bottom
    if is_reverse
        dists1, dists2 = dists_to_bottom, dists_to_top
        mol_order = reverse(mol_order)
    end

    min_n_iters = min(length(mol_order), round(Int, 10 / (1 - passed_frac_threshold)))
    border = mol_order[1]
    n_passed = 0.0
    for (it, mi) in enumerate(mol_order)
        n_passed += (dists1[mi] < dists2[mi])
        if (n_passed / it) >= passed_frac_threshold
            border = mi
        elseif it >= min_n_iters
            break
        end
    end

    return border
end

function estimate_noise_level(pos_data::Matrix{Float64}, mean_dists::Vector{Float64}; min_noise_frac::Float64=0.1, min_real_frac::Float64=0.2,
        max_iters::Int=100, quantile_change_threshold::Float64=0.01)#::Tuple{Float64, Float64}
    current_qs = quantile(mean_dists, (min_real_frac, 1.0 - min_noise_frac))
    for it in 1:max_iters
        top_q, bottom_q = current_qs
        top_mols = findall(mean_dists .<= top_q);
        bottom_mols = findall(mean_dists .>= bottom_q);
        unknown_ids = findall((mean_dists .> top_q) .& (mean_dists .< bottom_q));

        dists_to_top = getindex.(knn(KDTree(pos_data[:, top_mols]), pos_data[:, unknown_ids], 1, true)[2], 1);
        dists_to_bottom = getindex.(knn(KDTree(pos_data[:, bottom_mols]), pos_data[:, unknown_ids], 1, true)[2], 1);
        unknown_dists = mean_dists[unknown_ids];
        unk_order = sortperm(unknown_dists)

        old_qs = current_qs
        current_qs = [unknown_dists[find_noise_border(dists_to_top, dists_to_bottom, unk_order, is_reverse=r)] for r in [false, true]]
        if maximum((current_qs .- old_qs) ./ max.(current_qs, old_qs)) < quantile_change_threshold
            break
        end
    end

    return current_qs
end

"""
Parameters:
- min_noise_level: minimal fraction of molecules, which are hard assigned to the noise class
- max_noise_level: maximal fraction of molecules, which are hard assigned to the noise class
- min_real_level: minimal fraction of molecules, which are hard assigned to the real class
"""
function append_confidence!(df_spatial::DataFrame, segmentation_mask::Union{BitArray{2}, Nothing}=nothing; nn_id::Int,
        min_noise_frac::Float64=0.1, min_real_frac::Float64=0.1, max_noise_frac::Float64=0.6,
        seg_quantile_left::Float64=0.9, seg_quantile_mid::Float64=0.99)
    pos_data = position_data(df_spatial);
    mean_dists = getindex.(knn(KDTree(pos_data), pos_data, nn_id + 1, true)[2], nn_id + 1)
    border_left, border_right = estimate_noise_level(pos_data, mean_dists; min_noise_frac=min_noise_frac, min_real_frac=min_real_frac)
    border_right = max(border_right, 1 - max_noise_frac)

    if segmentation_mask !== nothing
        mol_in_nuclei_mask = segmentation_mask[CartesianIndex.(round.(Int, pos_data[2,:]), round.(Int, pos_data[1,:]))]
        sb_left, sb_mid = quantile(mean_dists[mol_in_nuclei_mask], [seg_quantile_left, seg_quantile_mid])
        border_left = max(sb_left, border_left) # TODO: maybe estimate with DAPI will be more robust, so should remove min/max
        border_right = min(sb_left + 2 * (sb_mid - sb_left), border_right)
    end

    df_spatial[!,:confidence] = interpolate_linear.(mean_dists, border_left, border_right);
end