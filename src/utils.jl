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

function split(array::Array{Int, 1}, factor::Array{Int, 1}; max_factor::Union{Int, Nothing}=nothing)
    @assert length(array) == length(factor)
    if max_factor === nothing
        max_factor = maximum(factor)
    end

    splitted = [Int[] for i in 1:max_factor]
    for i in 1:length(array)
        push!(splitted[factor[i]], array[i])
    end

    return splitted
end

split(array::UnitRange{Int64}, factor::Array{Int64,1}) = split(collect(array), factor)