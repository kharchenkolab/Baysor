using DataFrames
using KernelDensity

estimate_density_kde(coords::Array{Float64, 2}, points::Array{Float64, 2}, bandwidth::T where T <: Real)::Vector{Float64} =
    InterpKDE(kde((coords[1,:], coords[2,:]), bandwidth=(Float64(bandwidth), Float64(bandwidth)))).itp.(points[1,:], points[2,:])

function count_array(values::Union{Array{Int, 1}, SubArray{Int,1}}; max_value::Union{Int, Nothing}=nothing)
    if max_value === nothing
        max_value = maximum(values)
    end

    counts = zeros(Int, max_value)
    for v in values
        counts[v] += 1
    end

    return counts
end

function count_array!(counts::T1 where T1 <: AbstractArray{Int, 1}, values::T2 where T2 <: AbstractArray{Int, 1})
    counts .= 0
    for v in values
        counts[v] += 1
    end

    return counts
end

function prob_array(values::Union{Array{Int, 1}, SubArray{Int,1}}; max_value::Union{Int, Nothing}=nothing, smooth::Float64=0.0)
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

function prob_array!(counts::Union{Array{Float64, 1}, SubArray{Float64,1}}, values::Array{Int, 1}; smooth::Float64=0.0)
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

trim_mean(x::Array{T, 1} where T <: Real; prop::Float64=0.2) = mean(trim(x; prop=prop))

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

split(array::UnitRange{Int64}, factor::Array{Int64,1}; kwargs...) = split(collect(array), factor; kwargs...)

split_ids(factor::Array{Int, 1}; kwargs...) = split(1:length(factor), factor; kwargs...)

split(df::DataFrame, factor::Symbol; kwargs...) = split(df, Array(df[!, factor]); kwargs...)
split(df::DataFrame, factor::Array{Int, 1}; kwargs...) = [df[ids, :] for ids in split(1:size(df, 1), factor; kwargs...)]

function interpolate_linear(x::T, x_start::T, x_end::T; y_start::T=1.0, y_end::T=0.0)::Float64 where T<:Real
    if x < x_start
        return y_start
    elseif x > x_end
        return y_end
    end

    return y_start + (x - x_start) / (x_end - x_start) * (y_end - y_start)
end

function is_point_in_polygon(point::Union{Vector{T1}, Tuple{T1, T1}} where T1 <: Real, poly::Array{Vector{T2}, 1} where T2 <: Real,
        borders::Union{Array{Vector{T3},1}, Nothing} where T3 = nothing)
    if borders !== nothing && (borders[1][1] > point[1] || borders[1][2] > point[2] || borders[2][1] < point[1] || borders[2][2] < point[2])
        return false
    end

    j = length(poly)
    c = false
    for i in 1:length(poly)
        if ((poly[i][2] > point[2]) != (poly[j][2] > point[2])) &&
            (point[1] < poly[i][1] + (poly[j][1] - poly[i][1]) * (point[2] - poly[i][2]) / (poly[j][2] - poly[i][2]))
            c = !c
        end

        j = i
    end

    return c
end

"""Golden section search
to find the minimum of f on [a,b]
opt_func: a strictly unimodal function on [a,b]
"""
function linsearch_gs(opt_func::Function, a::T, b::T; tol=1e-3) where T<: Real
    gr = (sqrt(5) + 1) / 2

    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(c - d) > tol
        if opt_func(c) < opt_func(d)
            b = d
        else
            a = c
        end
        # We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
    end

    return (b + a) / 2, opt_func((b + a) / 2)
end

function trace_values_along_line(arr::Matrix{T}, start_x::Int, start_y::Int, end_x::Int, end_y::Int)::Vector{T} where T <: Real
    @assert all([end_y, end_x] .<= size(arr))
    @assert all([start_y, start_x] .<= size(arr))
    @assert all([start_y, start_x] .> 0)
    @assert all([end_y, end_x] .> 0)

    a = (end_y - start_y) / (end_x - start_x)
    b = end_y - a * end_x

    dx = sign(end_x - start_x)
    dy = sign(end_y - start_y)

    x, y = start_x, start_y

    vals = [arr[start_y, start_x]]
    while (x != end_x) | (y != end_y)
        if (abs((a * x + b) - (y + dy)) < abs((a * (x + dx) + b) - y)) || dx == 0
            y += dy
        else
            x += dx
        end
        push!(vals, arr[y, x])
    end
    return vals
end

@inline @fastmath function fmax(v1::Float64, v2::Float64)
    v1 > v2 ? v1 : v2
end

@inline @fastmath function fmin(v1::Float64, v2::Float64)
    v1 < v2 ? v1 : v2
end

### Statistics

@inline function fsample(w::Vector{Float64})::Int
    n = length(w)
    if n == 0
        error("Empty vector for sampling")
    end

    t = rand(Random.GLOBAL_RNG) * sum(w)
    i = 1
    cw = w[1]
    while cw < t && i < length(w)
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

@inline @inbounds fsample(arr::Vector{Int}, w::Vector{Float64})::Int = arr[fsample(w)]

"""
It works only for large samples with more or less uniform weights
"""
function wmedian(values::Vector{T} where T <: Real, weights::Vector{Float64}; ord::Vector{Int}=sortperm(values))
    w_avg = sum(weights) / 2
    w_cur = 0.0
    for i in ord
        w_cur += weights[i]
        if w_cur >= w_avg
            if i == 1
                return values[1]
            end

            return values[i]
        end
    end
end

wmad(values::Vector{T} where T <: Real, weights::Vector{Float64}; ord::Vector{Int}=sortperm(values)) =
    wmad(values, weights, wmedian(values, weights, ord=ord))

"""
It works only for large samples with more or less uniform weights
"""
wmad(values::Vector{T} where T <: Real, weights::Vector{Float64}, m::Float64) =
    1.4826 * wmedian(abs.(values .- m), weights)