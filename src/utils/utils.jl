using DataFrames
using KernelDensity

estimate_density_kde(coords::Array{Float64, 2}, points::Array{Float64, 2}, bandwidth::T where T <: Real)::Vector{Float64} =
    InterpKDE(kde((coords[1,:], coords[2,:]), bandwidth=(Float64(bandwidth), Float64(bandwidth)))).itp.(points[1,:], points[2,:])

function val_range(arr::AT where AT <: AbstractVector{<:Real})
    if length(arr) == 0
        return (nothing, nothing)
    end

    min_val, max_val = arr[1], arr[1]
    for v in arr
        min_val = fmin(min_val, v)
        max_val = fmax(max_val, v)
    end

    return min_val, max_val
end

count_array(values::VT where VT<: AbstractVector{<:Integer}, args...; max_value::Union{<:Integer, Nothing}=nothing, kwargs...) =
    count_array!(zeros(Int, something(max_value, maximum(values))), values, args...; erase_counts=false, kwargs...)

# TODO: check all usages and remove outer processing of 0 values
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

function split(vector::T where T <: AbstractVector; n_parts::Int)
    offset = ceil(Int, length(vector) / n_parts)
    return [vector[n:min(n + offset - 1, length(vector))] for n in 1:offset:length(vector)]
end

trim_mean(x::Array{T, 1} where T <: Real; prop::Float64=0.2) = mean(trim(x; prop=prop))

function split(array::T where T <: AbstractVector{TV}, factor::T2 where T2 <: AbstractVector{<:Integer}; max_factor::Union{Int, Nothing}=nothing)::Array{Vector{TV}, 1} where TV
    @assert length(array) == length(factor)
    if max_factor === nothing
        max_factor = maximum(factor)
    end

    splitted = [TV[] for i in 1:max_factor]
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

Returns: tuple with the optimal parameter and the optimal function value
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

@inline @fastmath function fmax(v1::T, v2::T) where T <: Real
    v1 > v2 ? v1 : v2
end

@inline @fastmath function fmin(v1::T, v2::T) where T <: Real
    v1 < v2 ? v1 : v2
end

function estimate_difference_l0(m1::Matrix{Float64}, m2::Matrix{Float64}; col_weights::Union{Nothing, Vector{Float64}}=nothing, change_threshold::Float64=1e-7)::Tuple{Float64, Float64}
    max_diff = 0.0
    n_changed = 0
    if !all(size(m1) .== size(m2))
        error("Matrices must be of the same size")
    end

    @inbounds for ci in 1:size(m1, 2)
        c_max = 0.0
        for ri in 1:size(m1, 1)
            c_max = fmax(abs(m1[ri, ci] - m2[ri, ci]), c_max)
        end

        if col_weights !== nothing
            c_max *= col_weights[ci]
        end

        if c_max > change_threshold
            n_changed += 1
        end

        max_diff = fmax(c_max, max_diff)
    end

    return max_diff, n_changed / size(m1, 2)
end

function pairwise_jaccard(values::Array{Vector{Int}, 1}, min_dist::Float64=0.0001)::Matrix{Float64}
    dist_mat = zeros(length(values), length(values))
    for i1 in 1:length(values)
        s1 = values[i1]
        for i2 in (i1+1):length(values)
            s2 = values[i2]
            inter_len = 0
            for v in s1
                if v in s2
                    inter_len += 1
                end
            end
            dist_mat[i1, i2] = dist_mat[i2, i1] = fmax(1.0 - length(inter_len) / (length(s1) + length(s2) - inter_len), min_dist)
        end
    end

    return dist_mat
end

function sparse_counts(vec1::Vector{<:Integer}, vec2::Vector{<:Integer})
    counts = Dict{Tuple{Int, Int}, Int}()
    for tp in zip(vec1, vec2)
        counts[tp] = get(counts, tp, 0) + 1
    end
    i1 = getindex.(keys(counts), 1)
    i2 = getindex.(keys(counts), 2)
    v = collect(values(counts))
    return sparse(i1, i2, v)
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
function wmedian(values::Vector{T} where T <: Real, weights::Vector{Float64}; ord::Vector{Int}=sortperm(values))::Float64
    w_avg = sum(weights) / 2
    w_cur = 0.0
    for i in ord
        w_cur += weights[i]
        if w_cur >= w_avg
            return values[i]
        end
    end

    if any(isnan.(weights))
        error("NaNs are presented in weights for wmedian")
    end
end

wmad(values::Vector{T} where T <: Real, weights::Vector{Float64}; ord::Vector{Int}=sortperm(values))::Float64 =
    wmad(values, weights, wmedian(values, weights, ord=ord))

"""
It works only for large samples with more or less uniform weights
"""
wmad(values::Vector{T} where T <: Real, weights::Vector{Float64}, m::Float64)::Float64 =
    1.4826 * wmedian(abs.(values .- m), weights)

function wmean(values::Vector{T} where T <: Real, weights::T1 where T1 <: AbstractVector{Float64}; non_zero_ids::T2 where T2 <: AbstractArray{Int} = 1:length(values))
    s, ws = 0.0, 0.0
    for i in non_zero_ids
        s += values[i] * weights[i]
        ws += weights[i]
    end

    return s / ws
end

function wmean_std(values::Vector{T0} where T0 <: Real, weights::T1 where T1 <: AbstractVector{Float64}; non_zero_ids::T2 where T2 <: AbstractArray{Int} = 1:length(values))
    m = wmean(values, weights; non_zero_ids=non_zero_ids)
    s, ws = 0.0, 0.0
    for i in non_zero_ids
        dv = (values[i] - m)
        s += dv * dv * weights[i]
        ws += weights[i]
    end

    return m, sqrt(s / ws)
end