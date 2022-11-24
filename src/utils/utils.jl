count_array(values::VT where VT<: AbstractVector{<:Integer}, args...; max_value::Union{<:Integer, Nothing}=nothing, kwargs...) =
    count_array!(zeros(Int, max_value !== nothing ? max_value : maximum(values)), values, args...; erase_counts=false, kwargs...)

function count_array!(counts::VT1 where VT1 <: AbstractVector{<:Integer}, values::VT2 where VT2 <: AbstractVector{<:Integer}; drop_zero::Bool=false, erase_counts::Bool=true)
    if erase_counts
        counts .= 0
    end

    has_zero = false
    for v in values
        if v == 0
            has_zero = true
            continue
        end

        counts[v] += 1
    end

    if !drop_zero && has_zero
        @warn "Array has zero values. It was ignored."
    end

    return counts
end

function count_array!(counts::VT1 where VT1 <: AbstractVector{RT}, values::VT2 where VT2 <: AbstractVector{<:Integer}, weights::VT3 where VT3 <: AbstractVector{RT};
        drop_zero::Bool=false, erase_counts::Bool=true) where RT<:Real
    if erase_counts
        counts .= 0
    end

    has_zero = false
    for i in eachindex(values)
        v = values[i]
        if v == 0
            has_zero = true
            continue
        end

        counts[v] += weights[i]
    end

    if !drop_zero && has_zero
        @warn "Array has zero values. It was ignored."
    end

    return counts
end

@inline @fastmath function fmax(v1::T, v2::T) where T <: Real
    v1 > v2 ? v1 : v2
end

@inline @fastmath function fmin(v1::T, v2::T) where T <: Real
    v1 < v2 ? v1 : v2
end

function default_param_value(
        param::Symbol, min_molecules_per_cell::Union{Int, Nothing};
        n_molecules::Union{Int, Nothing}=nothing, n_genes::Union{Int, Nothing}=nothing
    )
    if min_molecules_per_cell === nothing
        error("Either `$param` or `min_molecules_per_cell` must be provided")
    end

    min_molecules_per_cell = max(min_molecules_per_cell, 3)

    if param == :min_molecules_per_segment
        return max(round(Int, min_molecules_per_cell / 4), 2)
    end

    if param == :confidence_nn_id
        return max(div(min_molecules_per_cell, 2) + 1, 5)
    end

    if param == :composition_neighborhood
        if n_genes === nothing
            return max(min_molecules_per_cell, 3)
        end

        return max(div(n_genes, 10), min_molecules_per_cell, 3)
    end

    if param == :n_gene_pcs
        (n_genes !== nothing) || error("Either `$param` or `n_genes` must be provided")
        return min(max(div(n_genes, 3), 30), 100, n_genes)
    end

    if param == :n_cells_init
        (n_molecules !== nothing) || error("Either `$param` or `n_molecules` must be provided")
        return div(n_molecules, min_molecules_per_cell) * 2
    end
end