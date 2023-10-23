using SparseArrays
using Base.Threads
using NearestNeighbors: KDTree, knn

function knn_parallel(nn_tree::KDTree, x::AbstractMatrix{<:Real}, nn_interpolate::Int; sorted::Bool=false)
    indices = Vector{Vector{Int}}(undef, size(x, 2))
    distances = Vector{Vector{eltype(eltype(nn_tree.data))}}(undef, size(x, 2))

    @threads for i in axes(x, 2)
        indices[i], distances[i] = knn(nn_tree, x[:, i], nn_interpolate, sorted)
    end

    return indices, distances
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

struct PseudoWeight end;
Base.getindex(::PseudoWeight, i::Int) = 1
is_provided(::PseudoWeight) = false
is_provided(::Vector) = true

count_array_sparse(values::AbstractVector{Union{Missing, Int}}, args...; kwargs...) =
    count_array_sparse(collect(skipmissing(values)), args...; kwargs...)

count_array_sparse(T::DataType, values::AbstractVector{Union{Missing, Int}}, args...; kwargs...) =
    count_array_sparse(T, collect(skipmissing(values)), args...; kwargs...)

count_array_sparse(values::AbstractVector{Int}, ::Nothing=nothing; kwargs...) = count_array_sparse(Int, values; kwargs...)
count_array_sparse(values::AbstractVector{Int}, weights::AbstractVector{Float64}; kwargs...) =
    count_array_sparse(Float64, values, weights; kwargs...)

function count_array_sparse(
        T::DataType, values::AbstractVector{Int}, weights::Union{AbstractVector{Float64}, PseudoWeight}=PseudoWeight();
        total::Int=0, min_val::Float64=1e-5, normalize::Bool=false
    )
    !isempty(values) || return spzeros(T, total)

    if total <= 0
        total = maximum(values)
    end

    indices = Int[]
    counts = T[]

    last_id = values[1]
    cnt = 0
    id = 0
    for i in sortperm(values)
        id = values[i]
        if id != last_id
            if cnt > min_val
                push!(counts, cnt)
                push!(indices, last_id)
            end
            cnt = 0
            last_id = id
        end
        cnt += weights[i]
    end

    if cnt > min_val
        push!(counts, cnt)
        push!(indices, id)
    end

    if normalize
        counts = counts ./ sum(counts)
    end

    return SparseVector(total, indices, counts)
end

split(df::DataFrame, factor::Symbol; kwargs...) = split(df, Array(df[!, factor]); kwargs...)
split(df::DataFrame, factor::Vector{Int}; kwargs...) = [df[ids, :] for ids in split(1:size(df, 1), factor; kwargs...)]

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

### Statistics

@inline function fsample(w::AbstractVector{Float64})::Int
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

@inline @inbounds fsample(arr::Vector{Int}, w::AbstractVector{Float64})::Int = arr[fsample(w)]

function wmean(
        values::AbstractVector{<:Real}, weights::T where T <: AbstractVector{<:Real}
    )
    s, ws = 0.0, 0.0
    for i in eachindex(values)
        s += values[i] * weights[i]
        ws += weights[i]
    end

    return s / ws
end

function wmean_std(values::AbstractVector{T}, weights::AbstractVector{T}) where T<:Real
    # Somehow, for views to work fast, we need to estimate everything in a single pass
    # So I use an incremental algorithm here.

    # Initialize variables
    w = zero(T)
    w_sq = zero(T)
    w_x = zero(T)
    w_x_sq = zero(T)

    # Calculate weighted sums
    for i in eachindex(values)
        w_i = weights[i]
        x_i = values[i]

        w += w_i
        w_sq += w_i^2
        w_x += w_i*x_i
        w_x_sq += w_i*x_i^2
    end

    # Calculate weighted mean and variance
    m = w_x/w
    s2 = (w_x_sq/w) - m^2

    # Calculate weighted standard deviation
    s = sqrt(s2)

    return m,s
end

function sphstack(X::Vector{SparseVector{Tv, Ti}}) where {Tv, Ti}
    # Implementation copied from SparseArrays.jl
    # https://github.com/JuliaSparse/SparseArrays.jl/blob/f8f0f409972b627e05a5e89f58f189f94e3a5666/src/sparsevector.jl#L1114
    # See https://github.com/JuliaLang/julia/issues/50623 for the reason
    n = length(X)
    m = length(X[1])
    tnnz = nnz(X[1])
    for j = 2:n
        length(X[j]) == m || throw(DimensionMismatch("Inconsistent column lengths."))
        tnnz += nnz(X[j])
    end

    # construction
    colptr = Vector{Ti}(undef, n+1)
    nzrow = Vector{Ti}(undef, tnnz)
    nzval = Vector{Tv}(undef, tnnz)
    roff = 1
    @inbounds for j = 1:n
        xj = X[j]
        xnzind = SparseArrays.nonzeroinds(xj)
        xnzval = SparseArrays.nonzeros(xj)
        colptr[j] = roff
        copyto!(nzrow, roff, xnzind)
        copyto!(nzval, roff, xnzval)
        roff += length(xnzind)
    end
    colptr[n+1] = roff
    return SparseMatrixCSC{Tv,Ti}(m, n, colptr, nzrow, nzval)
end