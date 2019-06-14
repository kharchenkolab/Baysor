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

function split(df::DataFrame, factor::Array{Int, 1})
    res = Array{DataFrame, 1}(undef, maximum(factor))
    for i in unique(factor)
        res[i] = df[factor .== i, :]
    end

    return res
end

function convex_hull(a::Array{Array{T,1},1} where T<:Real)
    cw(a, b, c) = (a[1]*(b[2]-c[2])+b[1]*(c[2]-a[2])+c[1]*(a[2]-b[2]) < 0);
    ccw(a, b, c) = (a[1]*(b[2]-c[2])+b[1]*(c[2]-a[2])+c[1]*(a[2]-b[2]) > 0);

    if length(a) == 1
        return a
    end

    a = sort(a, lt = (a,b) -> (a[1] < b[1] || a[1] == b[1] && a[2] < b[2]))
    p1 = a[1];
    p2 = a[end];
    up = [p1]
    down = [p1]

    for i in 2:length(a)
        if (i==length(a) || cw(p1, a[i], p2))
            while (length(up) >= 2 && !cw(up[end-1], up[end], a[i]))
                up = up[1:end-1]
            end

            push!(up, a[i])
        end
        if (i==length(a) || ccw(p1, a[i], p2))
            while (length(down) >= 2 && !ccw(down[end-1], down[end], a[i]))
                down = down[1:end-1]
            end

            push!(down, a[i])
        end
    end

    return hcat(vcat(up, reverse(down[1:end-1]))...);
end

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene)
    df_spatial = CSV.read(data_path);

    if gene_col === nothing
        df_spatial = df_spatial[[x_col, y_col]]
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y);
    else
        df_spatial = df_spatial[[x_col, y_col, gene_col]]
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y, gene_col => :gene);
    end

    return df_spatial
end