using DataFrames

function count_array(values::Array{Int, 1}; max_value::Union{Int, Nothing}=nothing)
    if max_value === nothing
        max_value = maximum(values)
    end

    counts = zeros(Int, max_value)
    for v in values
        counts[v] += 1
    end

    return counts
end

function count_array!(counts::Array{Int, 1}, values::Array{Int, 1})
    counts .= 0
    for v in values
        counts[v] += 1
    end

    return counts
end

function prob_array(values::Array{Int, 1}; max_value::Union{Int, Nothing}=nothing, smooth::Float64=0.0)
    if max_value === nothing
        max_value = maximum(values)
    end

    sum_value = length(values) + max_value * smooth
    counts = fill(smooth / sum_value, max_value)
    for v in values
        counts[v] += 1.0 / sum_value
    end

    return counts
end

function prob_array!(counts::Array{Float64, 1}, values::Array{Int, 1}; smooth::Float64=0.0)
    sum_value = length(values) + length(counts) * smooth
    counts .= smooth / sum_value
    for v in values
        counts[v] += 1.0 / sum_value
    end

    return counts
end

function split(vector::AbstractArray{T, 1} where T; n_parts::Int)
    offset = ceil(Int, length(vector) / n_parts)
    return [vector[n:min(n + offset - 1, length(vector))] for n in 1:offset:length(vector)]
end

trim_mean(x::Array{Float64, 1}; prop::Float64=0.2) = mean(trim(x; prop=prop))

function split(array::Array{T, 1}, factor::Array{Int, 1}; max_factor::Union{Int, Nothing}=nothing)::Array{Array{T, 1}, 1} where T
    @assert length(array) == length(factor)
    if max_factor === nothing
        max_factor = maximum(factor)
    end

    splitted = [T[] for i in 1:max_factor]
    for i in 1:length(array)
        push!(splitted[factor[i]], array[i])
    end

    return splitted
end

split(array::UnitRange{Int64}, factor::Array{Int64,1}) = split(collect(array), factor)

split(df::DataFrame, factor::Symbol) = split(df, Array(df[!, factor]))
split(df::DataFrame, factor::Array{Int, 1}) = [df[ids, :] for ids in split(1:size(df, 1), factor)]

function interpolate_linear(x::T, x_start::T, x_end::T; y_start::T=1.0, y_end::T=0.0)::Float64 where T<:Real
    if x < x_start
        return y_start
    elseif x > x_end
        return y_end
    end

    return y_start + (x - x_start) / (x_end - x_start) * (y_end - y_start)
end