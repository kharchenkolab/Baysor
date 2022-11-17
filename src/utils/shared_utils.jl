# simple functions imported from different sub-modules

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
    for i in 1:length(values)
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